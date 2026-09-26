import hashlib
import io
import importlib.util
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
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


class GctaSetupTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.destination = self.root / "tools" / "gcta" / "1.94.1" / "gcta64"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def create_package(self) -> Path:
        archive = self.root / "gcta.tar.bz2"
        executable = b"\x7fELFtest-gcta"
        with tarfile.open(archive, mode="w:bz2") as package:
            member = tarfile.TarInfo(finemap.GCTA_ARCHIVE_MEMBER)
            member.size = len(executable)
            package.addfile(member, io.BytesIO(executable))
        return archive

    def test_install_gcta_verifies_and_extracts_project_local_binary(self):
        archive = self.create_package()
        expected_digest = hashlib.sha256(archive.read_bytes()).hexdigest()

        def copy_package(url, filename):
            self.assertEqual(url, finemap.GCTA_PACKAGE_URL)
            Path(filename).write_bytes(archive.read_bytes())

        with (
            mock.patch.object(finemap.platform, "system", return_value="Linux"),
            mock.patch.object(finemap.platform, "machine", return_value="x86_64"),
            mock.patch.object(
                finemap.urllib.request, "urlretrieve", side_effect=copy_package
            ),
            mock.patch.object(finemap, "GCTA_PACKAGE_SHA256", expected_digest),
        ):
            installed = finemap.install_gcta(self.destination)

        self.assertEqual(installed, self.destination.resolve())
        self.assertEqual(installed.read_bytes(), b"\x7fELFtest-gcta")
        self.assertTrue(installed.stat().st_mode & 0o111)

    def test_resolve_gcta_binary_prefers_path_executable(self):
        with (
            mock.patch.object(
                finemap.shutil, "which", return_value="/usr/bin/gcta64"
            ),
            mock.patch.object(finemap, "install_gcta") as install_mock,
        ):
            executable = finemap.resolve_gcta_binary()

        self.assertEqual(executable, "/usr/bin/gcta64")
        install_mock.assert_not_called()

    def test_resolve_gcta_binary_installs_missing_default(self):
        with (
            mock.patch.object(finemap.shutil, "which", return_value=None),
            mock.patch.object(
                finemap, "install_gcta", return_value=self.destination
            ) as install_mock,
        ):
            executable = finemap.resolve_gcta_binary()

        self.assertEqual(executable, str(self.destination))
        install_mock.assert_called_once_with()

    def test_resolve_gcta_binary_rejects_missing_override(self):
        with (
            mock.patch.object(finemap.shutil, "which", return_value=None),
            self.assertRaisesRegex(FileNotFoundError, "requested GCTA binary"),
        ):
            finemap.resolve_gcta_binary("/missing/gcta64")


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
            input_file = (
                self.gwas_dir
                / f"lifespan_female_{chrom}.{finemap.PHENO_NAME}.glm.linear"
            )
            pd.DataFrame(rows).to_csv(input_file, sep="\t", index=False)

        result = finemap.merge_gwas_sumstats()

        self.assertEqual(len(result), len(finemap.FLY_CHROMS))
        self.assertEqual(result.columns.tolist(), finemap.GWAS_COLUMNS)
        self.assertEqual(set(result["CHR"]), set(finemap.FLY_CHROMS))
        self.assertFalse(result[finemap.NUMERIC_COLUMNS].isna().any(axis=None))

    def test_merge_gwas_sumstats_rejects_missing_columns(self):
        input_file = self.gwas_dir / f"lifespan_female_2L.{finemap.PHENO_NAME}.glm.linear"
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
        self.susie_results_dir = self.root / "output" / "susie_results"
        self.ld_source.parent.mkdir()
        self.path_patch = mock.patch.multiple(
            finemap,
            MERGED_QC_SOURCE=self.ld_source,
            LD_REF_BFILE=self.ld_target,
            COJO_INPUT_FILE=self.cojo_input,
            COJO_OUT_PREFIX=self.cojo_prefix,
            REGIONS_DIR=self.regions_dir,
            SUSIE_WORK_DIR=self.susie_work_dir,
            SUSIE_RESULTS_DIR=self.susie_results_dir,
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
                "A1": ["A", "C"],
                "A2": ["G", "T"],
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

        invalid_alleles = region.copy()
        invalid_alleles.loc[0, "A2"] = invalid_alleles.loc[0, "A1"]
        invalid_alleles.to_csv(region_file, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "identical A1 and A2 alleles"):
            finemap.load_region_summary(region_file)

    def test_write_region_snplist_preserves_validated_snp_order(self):
        region_file = self.root / "chr2L_pos100_snps.tsv"
        pd.DataFrame(
            {
                "SNP": ["2L_50", "2L_100", "2L_150"],
                "A1": ["A", "C", "G"],
                "A2": ["G", "T", "A"],
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
                "A1": ["A", "C"],
                "A2": ["G", "T"],
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

    def test_compute_region_ld_runs_plink_and_reuses_nonempty_matrix(self):
        bfile_prefix = self.susie_work_dir / "chr2L_pos100"
        bfile_prefix.parent.mkdir(parents=True)
        for extension in (".bed", ".bim", ".fam"):
            bfile_prefix.with_suffix(extension).write_text("genotype\n", encoding="utf-8")
        ld_file = bfile_prefix.with_suffix(".ld")

        def create_ld(command, check):
            self.assertTrue(check)
            ld_file.write_text("1.0 0.5\n0.5 1.0\n", encoding="utf-8")

        with mock.patch.object(finemap.shutil, "which", return_value="/opt/plink"), mock.patch.object(
            finemap.subprocess, "run", side_effect=create_ld
        ) as run_mock:
            result = finemap.compute_region_ld(bfile_prefix, plink_bin="/opt/plink")
            cached_result = finemap.compute_region_ld(bfile_prefix, plink_bin="/opt/plink")

        self.assertEqual(result, ld_file)
        self.assertEqual(cached_result, ld_file)
        self.assertEqual(run_mock.call_count, 1)
        command = run_mock.call_args.args[0]
        self.assertEqual(command[0], "/opt/plink")
        self.assertEqual(command[command.index("--bfile") + 1], str(bfile_prefix))
        self.assertEqual(command[command.index("--r") + 1], "square")

    def test_load_aligned_region_data_orders_variants_and_flips_swapped_effects(self):
        region_file = self.root / "chr2L_pos100_snps.tsv"
        pd.DataFrame(
            {
                "SNP": ["2L_50", "2L_100"],
                "A1": ["A", "T"],
                "A2": ["G", "C"],
                "b": [0.1, 0.2],
                "se": [0.01, 0.02],
                "p": [1e-4, 1e-6],
                "N": [100, 101],
            }
        ).to_csv(region_file, sep="\t", index=False)
        bfile_prefix = self.susie_work_dir / "chr2L_pos100"
        bfile_prefix.parent.mkdir(parents=True)
        bfile_prefix.with_suffix(".bim").write_text(
            "1 2L_100 0 100 C T\n1 2L_50 0 50 A G\n", encoding="utf-8"
        )
        ld_file = bfile_prefix.with_suffix(".ld")
        ld_file.write_text("1.0 0.25\n0.25 1.0\n", encoding="utf-8")

        aligned, ld = finemap.load_aligned_region_data(region_file, bfile_prefix, ld_file)

        self.assertEqual(aligned["SNP"].tolist(), ["2L_100", "2L_50"])
        self.assertEqual(aligned["A1"].tolist(), ["C", "A"])
        self.assertEqual(aligned["b"].tolist(), [-0.2, 0.1])
        self.assertTrue(np.array_equal(ld, np.array([[1.0, 0.25], [0.25, 1.0]])))

        incompatible = pd.read_csv(region_file, sep="\t")
        incompatible.loc[0, ["A1", "A2"]] = ["C", "G"]
        incompatible.to_csv(region_file, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "alleles are incompatible"):
            finemap.load_aligned_region_data(region_file, bfile_prefix, ld_file)

    def test_run_susie_rss_returns_valid_pips_and_credible_sets(self):
        aligned = pd.DataFrame(
            {
                "b": [0.1, 0.2, -0.1],
                "se": [0.01, 0.02, 0.03],
                "N": [100, 102, 104],
            }
        )
        ld = np.array(
            [
                [1.0, 0.2, 0.1],
                [0.2, 1.0, 0.3],
                [0.1, 0.3, 1.0],
            ]
        )
        fake_ro = mock.MagicMock()
        fake_ro.NULL = object()
        fake_ro.FloatVector.side_effect = lambda values: list(values)
        fake_ro.r = {"matrix": mock.MagicMock(return_value="R_MATRIX")}

        fit = mock.MagicMock()
        fit.rx2.return_value = [0.8, 0.15, 0.05]
        cs_list = mock.MagicMock()
        cs_list.names = ["L1"]
        cs_list.rx2.return_value = [1, 2]
        cs_result = mock.MagicMock()
        cs_result.rx2.return_value = cs_list
        fake_susie = mock.MagicMock()
        fake_susie.susie_rss.return_value = fit
        fake_susie.susie_get_cs.return_value = cs_result

        pip, credible_sets = finemap.run_susie_rss(
            aligned, ld, runtime=(fake_ro, fake_susie)
        )

        self.assertTrue(np.array_equal(pip, np.array([0.8, 0.15, 0.05])))
        self.assertTrue(np.allclose(credible_sets[:2], [1.0, 1.0]))
        self.assertTrue(np.isnan(credible_sets[2]))
        call_kwargs = fake_susie.susie_rss.call_args.kwargs
        self.assertEqual(call_kwargs["n"], 102)
        self.assertEqual(call_kwargs["L"], finemap.SUSIE_MAX_EFFECTS)
        self.assertEqual(call_kwargs["R"], "R_MATRIX")
        fake_susie.susie_get_cs.assert_called_once_with(
            fit,
            coverage=finemap.SUSIE_COVERAGE,
            Xcorr="R_MATRIX",
        )

    def test_finemap_region_writes_sorted_results_and_reuses_output(self):
        region_file = self.root / "chr2L_pos100_snps.tsv"
        region_file.touch()
        bfile_prefix = self.susie_work_dir / "chr2L_pos100"
        ld_file = bfile_prefix.with_suffix(".ld")
        aligned = pd.DataFrame(
            {
                "SNP": ["2L_50", "2L_100", "2L_150"],
                "b": [0.1, 0.2, -0.1],
                "se": [0.01, 0.02, 0.03],
                "N": [100, 101, 102],
            }
        )
        ld = np.eye(3)

        with mock.patch.object(
            finemap, "extract_region_bfile", return_value=bfile_prefix
        ) as extract_mock, mock.patch.object(
            finemap, "compute_region_ld", return_value=ld_file
        ) as ld_mock, mock.patch.object(
            finemap, "load_aligned_region_data", return_value=(aligned, ld)
        ) as load_mock, mock.patch.object(
            finemap,
            "run_susie_rss",
            return_value=(np.array([0.1, 0.8, 0.1]), np.array([np.nan, 1.0, np.nan])),
        ) as susie_mock:
            output_file = finemap.finemap_region(region_file, plink_bin="/opt/plink")
            cached_file = finemap.finemap_region(region_file, plink_bin="/opt/plink")

        self.assertEqual(output_file, self.susie_results_dir / "chr2L_pos100_susie.tsv")
        self.assertEqual(cached_file, output_file)
        saved = pd.read_csv(output_file, sep="\t")
        self.assertEqual(saved["SNP"].tolist(), ["2L_100", "2L_150", "2L_50"])
        self.assertEqual(saved["PIP"].tolist(), [0.8, 0.1, 0.1])
        self.assertEqual(saved.loc[0, "CS"], 1.0)
        extract_mock.assert_called_once_with(region_file, plink_bin="/opt/plink", force=False)
        ld_mock.assert_called_once_with(bfile_prefix, plink_bin="/opt/plink", force=False)
        load_mock.assert_called_once_with(region_file, bfile_prefix, ld_file)
        susie_mock.assert_called_once()

    def test_run_susie_finemapping_reuses_runtime_and_removes_stale_results(self):
        region_files = [
            self.regions_dir / "chr2L_pos100_snps.tsv",
            self.regions_dir / "chrX_pos900_snps.tsv",
        ]
        output_files = [
            self.susie_results_dir / "chr2L_pos100_susie.tsv",
            self.susie_results_dir / "chrX_pos900_susie.tsv",
        ]
        self.susie_results_dir.mkdir(parents=True)
        stale_file = self.susie_results_dir / "chr3L_pos1_susie.tsv"
        stale_file.write_text("stale\n", encoding="utf-8")
        runtime = object()

        with mock.patch.object(
            finemap, "find_region_files", return_value=region_files
        ), mock.patch.object(
            finemap, "load_susie_runtime", return_value=runtime
        ) as runtime_mock, mock.patch.object(
            finemap, "finemap_region", side_effect=output_files
        ) as finemap_mock:
            results = finemap.run_susie_finemapping(plink_bin="/opt/plink")

        self.assertEqual(results, output_files)
        runtime_mock.assert_called_once_with()
        self.assertEqual(finemap_mock.call_count, 2)
        for call, region_file in zip(finemap_mock.call_args_list, region_files):
            self.assertEqual(call.args[0], region_file)
            self.assertEqual(call.kwargs["plink_bin"], "/opt/plink")
            self.assertFalse(call.kwargs["force"])
            self.assertIs(call.kwargs["runtime"], runtime)
        self.assertFalse(stale_file.exists())

    def test_run_cojo_stage_preserves_stage_one_execution_order(self):
        gwas = pd.DataFrame({"SNP": ["2L_100"]})
        significant_snps = pd.DataFrame({"SNP": ["2L_100"]})
        signals = pd.DataFrame({"SNP": ["2L_100"]})
        events = []

        with mock.patch.multiple(
            finemap,
            merge_gwas_sumstats=mock.DEFAULT,
            filter_significant_snps=mock.DEFAULT,
            write_cojo_input=mock.DEFAULT,
            prepare_cojo_bfile=mock.DEFAULT,
            run_cojo=mock.DEFAULT,
            load_cojo_signals=mock.DEFAULT,
            extract_regions=mock.DEFAULT,
        ) as stage_mocks:
            stage_mocks["merge_gwas_sumstats"].side_effect = lambda: events.append("merge") or gwas
            stage_mocks["filter_significant_snps"].side_effect = (
                lambda frame: events.append("filter") or significant_snps
            )
            stage_mocks["write_cojo_input"].side_effect = lambda frame: events.append("write")
            stage_mocks["prepare_cojo_bfile"].side_effect = lambda: events.append("reference")
            stage_mocks["run_cojo"].side_effect = lambda **kwargs: events.append("cojo")
            stage_mocks["load_cojo_signals"].side_effect = lambda: events.append("load") or signals
            stage_mocks["extract_regions"].side_effect = lambda frame, selected: events.append(
                "regions"
            )

            finemap.run_cojo_stage(gcta_bin="/opt/gcta64", force=True)

        self.assertEqual(
            events,
            ["merge", "filter", "write", "reference", "cojo", "load", "regions"],
        )
        stage_mocks["filter_significant_snps"].assert_called_once_with(gwas)
        stage_mocks["write_cojo_input"].assert_called_once_with(significant_snps)
        stage_mocks["run_cojo"].assert_called_once_with(gcta_bin="/opt/gcta64", force=True)
        stage_mocks["extract_regions"].assert_called_once_with(gwas, signals)

    def test_main_dispatches_selected_pipeline_stages(self):
        with mock.patch.object(finemap, "run_cojo_stage") as cojo_mock, mock.patch.object(
            finemap, "run_susie_finemapping"
        ) as susie_mock:
            finemap.main(
                [
                    "--stage",
                    "all",
                    "--gcta-bin",
                    "/opt/gcta64",
                    "--plink-bin",
                    "/opt/plink",
                    "--force",
                ]
            )

        cojo_mock.assert_called_once_with(gcta_bin="/opt/gcta64", force=True)
        susie_mock.assert_called_once_with(plink_bin="/opt/plink", force=True)

        with mock.patch.object(finemap, "run_cojo_stage") as cojo_mock, mock.patch.object(
            finemap, "run_susie_finemapping"
        ) as susie_mock:
            finemap.main(["--stage", "susie"])

        cojo_mock.assert_not_called()
        susie_mock.assert_called_once_with(plink_bin=finemap.PLINK_BIN, force=False)


if __name__ == "__main__":
    unittest.main()
