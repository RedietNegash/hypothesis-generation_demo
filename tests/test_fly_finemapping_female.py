import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "notebooks" / "fly_finemapping_female.py"
MODULE_SPEC = importlib.util.spec_from_file_location("fly_finemapping_female", MODULE_PATH)
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


if __name__ == "__main__":
    unittest.main()
