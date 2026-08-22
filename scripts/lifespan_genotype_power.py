#!/usr/bin/env python3

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import nct, t


BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")

BFILE = BASE_DIR / "data" / "magma" / "merged_qc"
PHENO_FILE = BASE_DIR / "data" / "gwas" / "lifespan_female.pheno"
GWAS_FILE = BASE_DIR / "data" / "gwas" / "lifespan_female.assoc"

OUTPUT_DIR = BASE_DIR / "data" / "magma" / "power"

RESULTS_FILE = OUTPUT_DIR / "genotype_power_results.tsv"
SUMMARY_FILE = OUTPUT_DIR / "genotype_power_summary.tsv"

ALPHA = 0.05

HYPOTHETICAL_DIFFERENCES = [
    2.0,
    5.0,
    10.0,
    15.0,
]

MAX_GWAS_SNPS = 100

PHENO_COLUMN = "S18_1537_F"


def read_phenotype():
    pheno = pd.read_csv(
        PHENO_FILE,
        sep=r"\s+",
        engine="python"
    )

    required = {"FID", "IID", PHENO_COLUMN}

    missing = required - set(pheno.columns)

    if missing:
        raise ValueError(
            f"Missing phenotype columns: {sorted(missing)}"
        )

    pheno = pheno[
        ["FID", "IID", PHENO_COLUMN]
    ].copy()

    pheno[PHENO_COLUMN] = pd.to_numeric(
        pheno[PHENO_COLUMN],
        errors="coerce"
    )

    pheno = pheno.dropna(
        subset=[PHENO_COLUMN]
    )

    pheno["FID"] = pheno["FID"].astype(str)
    pheno["IID"] = pheno["IID"].astype(str)

    return pheno


def read_fam():
    fam_file = BFILE.with_suffix(".fam")

    fam = pd.read_csv(
        fam_file,
        sep=r"\s+",
        header=None,
        names=[
            "FID",
            "IID",
            "PAT",
            "MAT",
            "SEX",
            "PHENO"
        ],
        dtype={
            "FID": str,
            "IID": str
        }
    )

    return fam


def read_gwas():
    if not GWAS_FILE.exists():
        raise FileNotFoundError(
            f"GWAS file not found: {GWAS_FILE}"
        )

    gwas = pd.read_csv(
        GWAS_FILE,
        sep=r"\s+",
        engine="python"
    )

    rename_map = {}

    for column in gwas.columns:
        upper = column.upper()

        if upper in {"SNP", "MARKER", "ID"}:
            rename_map[column] = "SNP"

        elif upper in {"P", "PVAL", "PVALUE", "P_VALUE"}:
            rename_map[column] = "P"

    gwas = gwas.rename(
        columns=rename_map
    )

    if "SNP" not in gwas.columns:
        raise ValueError(
            "Could not find SNP column in GWAS file."
        )

    if "P" not in gwas.columns:
        raise ValueError(
            "Could not find P-value column in GWAS file."
        )

    gwas["P"] = pd.to_numeric(
        gwas["P"],
        errors="coerce"
    )

    gwas = gwas.dropna(
        subset=["P"]
    )

    gwas = gwas.sort_values(
        "P"
    )

    return gwas.head(
        MAX_GWAS_SNPS
    )


def read_bim():
    bim_file = BFILE.with_suffix(".bim")

    bim = pd.read_csv(
        bim_file,
        sep=r"\s+",
        header=None,
        names=[
            "CHR",
            "SNP",
            "CM",
            "BP",
            "A1",
            "A2"
        ],
        dtype={
            "CHR": str,
            "SNP": str,
            "A1": str,
            "A2": str
        }
    )

    return bim


def calculate_cohens_d(
    mean1,
    mean2,
    var1,
    var2,
    n1,
    n2
):
    pooled_variance = (
        ((n1 - 1) * var1)
        + ((n2 - 1) * var2)
    ) / (
        n1 + n2 - 2
    )

    pooled_sd = math.sqrt(
        pooled_variance
    )

    if pooled_sd == 0:
        return np.nan

    return abs(
        mean1 - mean2
    ) / pooled_sd


def calculate_power(
    d,
    n1,
    n2,
    alpha=0.05
):
    if (
        not np.isfinite(d)
        or d <= 0
        or n1 < 2
        or n2 < 2
    ):
        return np.nan

    df = n1 + n2 - 2

    ncp = (
        d
        * math.sqrt(
            (n1 * n2)
            / (n1 + n2)
        )
    )

    critical = t.ppf(
        1 - alpha / 2,
        df
    )

    upper = 1 - nct.cdf(
        critical,
        df,
        ncp
    )

    lower = nct.cdf(
        -critical,
        df,
        ncp
    )

    return upper + lower


def get_snp_genotypes(snp):
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        snp_file = tmp / "snp.txt"
        snp_file.write_text(
            f"{snp}\n"
        )

        prefix = tmp / "subset"

        subprocess.run(
            [
                "plink",
                "--bfile",
                str(BFILE),
                "--extract",
                str(snp_file),
                "--recode",
                "A",
                "--out",
                str(prefix)
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        raw_file = prefix.with_suffix(
            ".raw"
        )

        if not raw_file.exists():
            return None

        raw = pd.read_csv(
            raw_file,
            sep=r"\s+",
            engine="python"
        )

        dosage_columns = [
            column
            for column in raw.columns
            if column not in {
                "FID",
                "IID",
                "PAT",
                "MAT",
                "SEX",
                "PHENOTYPE"
            }
        ]

        if not dosage_columns:
            return None

        dosage = raw[
            dosage_columns[0]
        ]

        result = pd.DataFrame(
            {
                "FID": raw["FID"].astype(str),
                "IID": raw["IID"].astype(str),
                "GENOTYPE": pd.to_numeric(
                    dosage,
                    errors="coerce"
                )
            }
        )

        return result


def analyze_snp(
    snp,
    phenotype
):
    genotype = get_snp_genotypes(
        snp
    )

    if genotype is None:
        return None

    data = phenotype.merge(
        genotype,
        on=["FID", "IID"],
        how="inner"
    )

    data = data.dropna(
        subset=[
            PHENO_COLUMN,
            "GENOTYPE"
        ]
    )

    if data.empty:
        return None

    genotype_counts = (
        data["GENOTYPE"]
        .value_counts()
        .sort_index()
    )

    available = [
        genotype
        for genotype, count
        in genotype_counts.items()
        if count >= 2
    ]

    if len(available) < 2:
        return None

    if 0 in available and 2 in available:
        g1 = 0
        g2 = 2
    else:
        g1 = available[0]
        g2 = available[-1]

    group1 = data[
        data["GENOTYPE"] == g1
    ][PHENO_COLUMN]

    group2 = data[
        data["GENOTYPE"] == g2
    ][PHENO_COLUMN]

    n1 = len(group1)
    n2 = len(group2)

    if n1 < 2 or n2 < 2:
        return None

    mean1 = group1.mean()
    mean2 = group2.mean()

    var1 = group1.var(
        ddof=1
    )

    var2 = group2.var(
        ddof=1
    )

    d = calculate_cohens_d(
        mean1,
        mean2,
        var1,
        var2,
        n1,
        n2
    )

    power = calculate_power(
        d,
        n1,
        n2,
        ALPHA
    )

    pooled_sd = (
        abs(mean1 - mean2) / d
        if np.isfinite(d) and d != 0
        else np.nan
    )

    return {
        "SNP": snp,
        "GENOTYPE_1": g1,
        "GENOTYPE_2": g2,
        "N1": n1,
        "N2": n2,
        "MEAN_1": mean1,
        "MEAN_2": mean2,
        "VAR_1": var1,
        "VAR_2": var2,
        "POOLED_SD": pooled_sd,
        "MEAN_DIFFERENCE": abs(
            mean1 - mean2
        ),
        "COHENS_D": d,
        "POWER": power
    }


def calculate_hypothetical_power(
    phenotype
):
    n = len(phenotype)

    results = []

    phenotype_sd = phenotype[
        PHENO_COLUMN
    ].std(
        ddof=1
    )

    for difference in HYPOTHETICAL_DIFFERENCES:
        d = (
            difference
            / phenotype_sd
        )

        power = calculate_power(
            d,
            n // 2,
            n - (n // 2),
            ALPHA
        )

        results.append(
            {
                "EFFECT_DAYS": difference,
                "ASSUMED_N1": n // 2,
                "ASSUMED_N2": n - (n // 2),
                "PHENOTYPE_SD": phenotype_sd,
                "COHENS_D": d,
                "POWER": power
            }
        )

    return pd.DataFrame(
        results
    )


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    print(
        "Loading phenotype..."
    )

    phenotype = read_phenotype()

    print(
        f"Phenotype samples: {len(phenotype)}"
    )

    print(
        "Loading GWAS results..."
    )

    gwas = read_gwas()

    print(
        f"SNPs evaluated: {len(gwas)}"
    )

    rows = []

    for i, snp in enumerate(
        gwas["SNP"],
        start=1
    ):
        print(
            f"[{i}/{len(gwas)}] {snp}"
        )

        try:
            result = analyze_snp(
                str(snp),
                phenotype
            )

            if result is not None:
                result["GWAS_P"] = float(
                    gwas.loc[
                        gwas["SNP"] == snp,
                        "P"
                    ].iloc[0]
                )

                rows.append(
                    result
                )

        except Exception as exc:
            print(
                f"Skipping {snp}: {exc}"
            )

    observed = pd.DataFrame(
        rows
    )

    observed.to_csv(
        RESULTS_FILE,
        sep="\t",
        index=False
    )

    hypothetical = (
        calculate_hypothetical_power(
            phenotype
        )
    )

    hypothetical.to_csv(
        SUMMARY_FILE,
        sep="\t",
        index=False
    )

    print(
        "\nObserved SNP power:"
    )

    if observed.empty:
        print(
            "No SNPs could be analyzed."
        )
    else:
        print(
            observed[
                [
                    "SNP",
                    "GWAS_P",
                    "N1",
                    "N2",
                    "MEAN_DIFFERENCE",
                    "COHENS_D",
                    "POWER"
                ]
            ].to_string(
                index=False
            )
        )

    print(
        "\nHypothetical effect-size power:"
    )

    print(
        hypothetical.to_string(
            index=False
        )
    )

    print(
        f"\nResults: {RESULTS_FILE}"
    )

    print(
        f"Summary: {SUMMARY_FILE}"
    )


if __name__ == "__main__":
    main()