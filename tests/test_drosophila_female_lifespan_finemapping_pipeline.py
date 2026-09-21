import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "notebooks"
    / "drosophila_female_lifespan_finemapping_pipeline.py"
)
MODULE_SPEC = importlib.util.spec_from_file_location(
    "drosophila_female_lifespan_finemapping_pipeline", MODULE_PATH
)
finemap = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(finemap)


class InputPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.gwas_dir = self.root / "gwas"
        self.output_dir = self.root / "finemap"
        self.gwas_dir.mkdir()
        self.path_patch = mock.patch.multiple(
            finemap,
            GLM_DIR=self.gwas_dir,
            OUT_DIR=self.output_dir,
            SIG_SNP_FILE=self.output_dir / "significant.tsv",
            COJO_INPUT_FILE=self.output_dir / "cojo.txt",
        )
        self.path_patch.start()

    def tearDown(self):
        self.path_patch.stop()
        self.temporary_directory.cleanup()

    def gwas_row(self, chrom, position, snp):
        return {
            "#CHROM": chrom,
            "POS": position,
            "ID": snp,
            "A1": "A",
            "OMITTED": "G",
            "A1_FREQ": 0.25,
            "BETA": 0.2,
            "SE": 0.05,
            "P": 1e-6,
            "OBS_CT": 150,
        }

    def test_merge_gwas_sumstats_standardizes_and_cleans_rows(self):
        for index, chrom in enumerate(finemap.FLY_CHROMS, start=1):
            rows = [self.gwas_row(chrom, index * 100, f"{chrom}_{index * 100}")]
            if chrom == "2L":
                invalid_row = self.gwas_row(chrom, 200, "2L_200")
                invalid_row["SE"] = "invalid"
                rows.append(invalid_row)
            input_file = self.gwas_dir / f"lifespan_{chrom}.{finemap.PHENO_NAME}.glm.linear"
            pd.DataFrame(rows).to_csv(input_file, sep="\t", index=False)

        result = finemap.merge_gwas_sumstats()

        self.assertEqual(len(result), len(finemap.FLY_CHROMS))
        self.assertEqual(result.columns.tolist(), finemap.GWAS_COLUMNS)
        self.assertEqual(set(result["CHR"]), set(finemap.FLY_CHROMS))
        self.assertFalse(result[finemap.NUMERIC_COLUMNS].isna().any(axis=None))

    def test_merge_gwas_sumstats_rejects_missing_columns(self):
        input_file = self.gwas_dir / f"lifespan_2L.{finemap.PHENO_NAME}.glm.linear"
        pd.DataFrame([self.gwas_row("2L", 100, "2L_100")]).drop(columns="SE").to_csv(
            input_file, sep="\t", index=False
        )

        with self.assertRaisesRegex(ValueError, "missing columns: SE"):
            finemap.merge_gwas_sumstats()

    def test_filter_significant_snps_applies_all_eligibility_rules(self):
        gwas = pd.DataFrame(
            {
                "CHR": ["2L"] * 5,
                "POS": [100, 200, 300, 400, 500],
                "SNP": ["kept", "low_maf", "low_n", "bad_p", "bad_se"],
                "A1": ["A"] * 5,
                "A2": ["G"] * 5,
                "freq": [0.05, 0.049, 0.2, 0.2, 0.2],
                "b": [0.1] * 5,
                "se": [0.01, 0.01, 0.01, 0.01, 0.0],
                "p": [1e-5, 1e-6, 1e-6, -1.0, 1e-6],
                "N": [100, 100, 99, 100, 100],
            }
        )

        result = finemap.filter_significant_snps(gwas)

        self.assertEqual(result["SNP"].tolist(), ["kept"])
        saved = pd.read_csv(finemap.SIG_SNP_FILE, sep="\t")
        self.assertEqual(saved["SNP"].tolist(), ["kept"])

    def test_write_cojo_input_preserves_required_format(self):
        significant_snps = pd.DataFrame(
            {
                "SNP": ["2L_100", "2L_200"],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
                "freq": [0.2, 0.3],
                "b": [0.1, 0.2],
                "se": [0.01, 0.02],
                "p": [1e-6, 2e-6],
                "N": [100.0, 101.0],
            }
        )

        finemap.write_cojo_input(significant_snps)

        saved = pd.read_csv(finemap.COJO_INPUT_FILE, sep=r"\s+")
        self.assertEqual(saved.columns.tolist(), finemap.COJO_COLUMNS)
        self.assertEqual(saved["N"].tolist(), [100, 101])

        duplicates = pd.concat([significant_snps, significant_snps.iloc[[0]]])
        with self.assertRaisesRegex(ValueError, "duplicate SNP IDs"):
            finemap.write_cojo_input(duplicates)


class PipelineOutputTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.ld_source = self.root / "source" / "merged_qc"
        self.ld_target = self.root / "output" / "bfile" / "merged_qc_numeric"
        self.cojo_input = self.root / "output" / "cojo_input.txt"
        self.cojo_prefix = self.root / "output" / "cojo" / "female_lifespan_cojo"
        self.regions_dir = self.root / "output" / "regions"
        self.susie_work_dir = self.root / "output" / "susie"
        self.ld_source.parent.mkdir()
        self.path_patch = mock.patch.multiple(
            finemap,
            MERGED_QC_SOURCE=self.ld_source,
            LD_REF_BFILE=self.ld_target,
            COJO_INPUT_FILE=self.cojo_input,
            COJO_OUT_PREFIX=self.cojo_prefix,
            REGIONS_DIR=self.regions_dir,
            SUSIE_WORK_DIR=self.susie_work_dir,
        )
        self.path_patch.start()

    def tearDown(self):
        self.path_patch.stop()
        self.temporary_directory.cleanup()

    def create_ld_reference(self):
        self.ld_source.with_suffix(".bed").write_bytes(b"bed")
        self.ld_source.with_suffix(".fam").write_text("sample\n", encoding="utf-8")
        self.ld_source.with_suffix(".bim").write_text(
            "2L rs1 0 100 A G\nX rs2 0 200 C T\n", encoding="utf-8"
        )

    def test_prepare_cojo_bfile_validates_and_remaps_chromosomes(self):
        self.create_ld_reference()

        finemap.prepare_cojo_bfile()

        self.assertTrue(self.ld_target.with_suffix(".bed").is_symlink())
        self.assertTrue(self.ld_target.with_suffix(".fam").is_symlink())
        remapped_rows = self.ld_target.with_suffix(".bim").read_text(encoding="utf-8").splitlines()
        self.assertEqual(remapped_rows[0].split("\t")[0], "1")
        self.assertEqual(remapped_rows[1].split("\t")[0], "23")

    def test_run_cojo_builds_expected_command_and_reuses_result(self):
        self.cojo_input.parent.mkdir(parents=True)
        self.cojo_input.write_text("input\n", encoding="utf-8")
        self.ld_target.parent.mkdir(parents=True)
        for extension in (".bed", ".bim", ".fam"):
            self.ld_target.with_suffix(extension).write_text("input\n", encoding="utf-8")

        result_file = self.cojo_prefix.with_suffix(".jma.cojo")
        result_file.parent.mkdir(parents=True)
        result_file.write_text("stale result\n", encoding="utf-8")

        def create_result(command, check):
            self.assertTrue(check)
            self.assertFalse(result_file.exists())
            result_file.write_text("fresh result\n", encoding="utf-8")

        with mock.patch.object(finemap.shutil, "which", return_value="/opt/gcta64"), mock.patch.object(
            finemap.subprocess, "run", side_effect=create_result
        ) as run_mock:
            finemap.run_cojo(gcta_bin="/opt/gcta64", force=True)
            finemap.run_cojo(gcta_bin="/opt/gcta64")

        self.assertEqual(run_mock.call_count, 1)
        self.assertEqual(result_file.read_text(encoding="utf-8"), "fresh result\n")
        command = run_mock.call_args.args[0]
        self.assertEqual(command[0], "/opt/gcta64")
        self.assertEqual(command[command.index("--bfile") + 1], str(self.ld_target))
        self.assertEqual(command[command.index("--cojo-file") + 1], str(self.cojo_input))
        self.assertEqual(command[command.index("--cojo-p") + 1], str(finemap.SIG_P_THRESHOLD))

    def test_load_cojo_signals_validates_and_sorts_results(self):
        result_file = self.cojo_prefix.with_suffix(".jma.cojo")
        result_file.parent.mkdir(parents=True)
        result_file.write_text(
            "Chr SNP bp b p bJ pJ\n"
            "23 rsX 900 0.2 1e-5 0.1 2e-5\n"
            "1 rs2L 100 0.3 1e-6 0.2 1e-6\n",
            encoding="utf-8",
        )

        signals = finemap.load_cojo_signals()

        self.assertEqual(signals["SNP"].tolist(), ["rs2L", "rsX"])
        self.assertEqual(signals["Chr"].tolist(), [1, 23])
        self.assertEqual(signals["bp"].tolist(), [100, 900])

        result_file.write_text(
            "Chr SNP bp b p bJ pJ\n6 invalid 100 0.3 1e-6 0.2 1e-6\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "unsupported chromosomes: 6"):
            finemap.load_cojo_signals()

    def test_extract_regions_sorts_outputs_and_removes_stale_files(self):
        self.regions_dir.mkdir(parents=True)
        stale_file = self.regions_dir / "chr3L_pos1_snps.tsv"
        stale_file.write_text("stale\n", encoding="utf-8")
        gwas = pd.DataFrame(
            {
                "CHR": ["2L", "2L", "X"],
                "POS": [150, 50, 900],
                "SNP": ["2L_150", "2L_50", "X_900"],
            }
        )
        signals = pd.DataFrame(
            {"Chr": [1, 23], "SNP": ["lead_2L", "lead_X"], "bp": [100, 900]}
        )

        finemap.extract_regions(gwas, signals)

        first_region = pd.read_csv(self.regions_dir / "chr2L_pos100_snps.tsv", sep="\t")
        self.assertEqual(first_region["SNP"].tolist(), ["2L_50", "2L_150"])
        self.assertTrue((self.regions_dir / "chrX_pos900_snps.tsv").is_file())
        self.assertFalse(stale_file.exists())

    def test_find_region_files_validates_and_orders_stage_two_inputs(self):
        self.regions_dir.mkdir(parents=True)
        with self.assertRaisesRegex(FileNotFoundError, "No fine-mapping region files"):
            finemap.find_region_files()

        second = self.regions_dir / "chrX_pos900_snps.tsv"
        first = self.regions_dir / "chr2L_pos100_snps.tsv"
        second.write_text("SNP\nX_900\n", encoding="utf-8")
        first.write_text("SNP\n2L_100\n", encoding="utf-8")
        (self.regions_dir / "unrelated.tsv").write_text("ignored\n", encoding="utf-8")

        region_files = finemap.find_region_files()

        self.assertEqual(region_files, [first, second])
        self.assertEqual(
            [finemap.region_label(path) for path in region_files],
            ["chr2L_pos100", "chrX_pos900"],
        )

    def test_load_region_summary_validates_susie_statistics(self):
        region_file = self.root / "region.tsv"
        region = pd.DataFrame(
            {
                "SNP": ["2L_100", "2L_200"],
                "b": [0.1, -0.2],
                "se": [0.01, 0.02],
                "p": [1e-6, 2e-5],
                "N": [100.0, 101.0],
            }
        )
        region.to_csv(region_file, sep="\t", index=False)

        loaded = finemap.load_region_summary(region_file)

        self.assertEqual(loaded["SNP"].tolist(), ["2L_100", "2L_200"])
        self.assertEqual(loaded["N"].tolist(), [100, 101])

        duplicated = pd.concat([region, region.iloc[[0]]], ignore_index=True)
        duplicated.to_csv(region_file, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "duplicate SNP IDs"):
            finemap.load_region_summary(region_file)

    def test_write_region_snplist_preserves_validated_snp_order(self):
        region_file = self.root / "chr2L_pos100_snps.tsv"
        pd.DataFrame(
            {
                "SNP": ["2L_50", "2L_100", "2L_150"],
                "b": [0.1, 0.2, 0.3],
                "se": [0.01, 0.02, 0.03],
                "p": [1e-4, 1e-6, 1e-5],
                "N": [100, 101, 102],
            }
        ).to_csv(region_file, sep="\t", index=False)

        snplist_file = finemap.write_region_snplist(region_file)

        self.assertEqual(snplist_file, self.susie_work_dir / "chr2L_pos100.snplist")
        self.assertEqual(
            snplist_file.read_text(encoding="utf-8").splitlines(),
            ["2L_50", "2L_100", "2L_150"],
        )

    def test_extract_region_bfile_runs_plink_and_reuses_complete_output(self):
        region_file = self.root / "chr2L_pos100_snps.tsv"
        pd.DataFrame(
            {
                "SNP": ["2L_50", "2L_100"],
                "b": [0.1, 0.2],
                "se": [0.01, 0.02],
                "p": [1e-4, 1e-6],
                "N": [100, 101],
            }
        ).to_csv(region_file, sep="\t", index=False)
        self.ld_target.parent.mkdir(parents=True)
        for extension in (".bed", ".bim", ".fam"):
            self.ld_target.with_suffix(extension).write_text("reference\n", encoding="utf-8")

        output_prefix = self.susie_work_dir / "chr2L_pos100"

        def create_bfile(command, check):
            self.assertTrue(check)
            output_prefix.with_suffix(".bed").write_bytes(b"bed")
            output_prefix.with_suffix(".bim").write_text(
                "1 2L_50 0 50 A G\n1 2L_100 0 100 C T\n", encoding="utf-8"
            )
            output_prefix.with_suffix(".fam").write_text("sample\n", encoding="utf-8")

        with mock.patch.object(finemap.shutil, "which", return_value="/opt/plink"), mock.patch.object(
            finemap.subprocess, "run", side_effect=create_bfile
        ) as run_mock:
            result = finemap.extract_region_bfile(region_file, plink_bin="/opt/plink")
            cached_result = finemap.extract_region_bfile(region_file, plink_bin="/opt/plink")

        self.assertEqual(result, output_prefix)
        self.assertEqual(cached_result, output_prefix)
        self.assertEqual(run_mock.call_count, 1)
        command = run_mock.call_args.args[0]
        self.assertEqual(command[0], "/opt/plink")
        self.assertEqual(command[command.index("--bfile") + 1], str(self.ld_target))
        self.assertEqual(
            command[command.index("--extract") + 1],
            str(self.susie_work_dir / "chr2L_pos100.snplist"),
        )


if __name__ == "__main__":
    unittest.main()
