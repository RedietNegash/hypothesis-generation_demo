#!/usr/bin/env python3

import math
import subprocess
from pathlib import Path

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
GWAS_DIR = BASE_DIR / "data" / "gwas"
TMP_DIR = GWAS_DIR / "tmp"
PHENO_FILE = GWAS_DIR / "lifespan_female.pheno"
INPUT_PREFIX = TMP_DIR / "merged_qc"
KEEP_FILE = TMP_DIR / "female_lifespan.keep"
GENO_PREFIX = TMP_DIR / "female_lifespan_genotypes"
FREQ_PREFIX = TMP_DIR / "female_lifespan_freq"
POWER_FILE = TMP_DIR / "female_lifespan_power.tsv"

N_EXPECTED = 197
ALPHA = 2.28e-8
EFFECTS = [1.0, 2.0, 5.0, 7.5, 10.0]
MAFS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]

def run(cmd):
    subprocess.run(cmd, check=True)

def verify_files():
    required = [
        INPUT_PREFIX.with_suffix(".bed"),
        INPUT_PREFIX.with_suffix(".bim"),
        INPUT_PREFIX.with_suffix(".fam"),
        PHENO_FILE,
    ]

    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

def create_keep_file():
    with open(PHENO_FILE) as f, open(KEEP_FILE, "w") as out:
        next(f)
        for line in f:
            fields = line.split()
            if len(fields) >= 2:
                out.write(f"{fields[0]} {fields[1]}\n")

def count_lines(path):
    with open(path) as f:
        return sum(1 for _ in f)

def calculate_phenotype_stats():
    values = []

    with open(PHENO_FILE) as f:
        next(f)
        for line in f:
            fields = line.split()
            if len(fields) >= 3:
                values.append(float(fields[2]))

    n = len(values)

    if n < 2:
        raise ValueError("Not enough phenotype observations")

    mean = sum(values) / n
    variance = sum((x - mean) ** 2 for x in values) / (n - 1)
    sd = math.sqrt(variance)

    return n, mean, sd

def create_genotype_subset():
    run([
        "plink",
        "--bfile", str(INPUT_PREFIX),
        "--keep", str(KEEP_FILE),
        "--allow-extra-chr",
        "--make-bed",
        "--out", str(GENO_PREFIX),
    ])

def calculate_maf():
    run([
        "plink",
        "--bfile", str(GENO_PREFIX),
        "--allow-extra-chr",
        "--freq",
        "--out", str(FREQ_PREFIX),
    ])

def get_maf_summary():
    mafs = []

    with open(FREQ_PREFIX.with_suffix(".frq")) as f:
        next(f)
        for line in f:
            fields = line.split()
            if len(fields) >= 5:
                maf = float(fields[4])
                if maf > 0:
                    mafs.append(maf)

    mafs.sort()

    if not mafs:
        raise ValueError("No valid MAF values found")

    n = len(mafs)
    median = mafs[n // 2]

    return len(mafs), mafs[0], median, mafs[-1]

def normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def power_for_effect(n, sd, maf, effect_days):
    genotype_variance = 2.0 * maf * (1.0 - maf)
    beta = effect_days / 2.0
    se = sd / math.sqrt(n * genotype_variance)

    z = abs(beta / se)
    z_alpha = 5.612001243305505

    power = (
        1.0 - normal_cdf(z_alpha - z)
        + normal_cdf(-z_alpha - z)
    )

    return min(max(power, 0.0), 1.0)

def write_power_table(n, sd):
    with open(POWER_FILE, "w") as out:
        out.write(
            "N\tMAF\teffect_days\tsd_lifespan\talpha\tpower\n"
        )

        for maf in MAFS:
            for effect in EFFECTS:
                power = power_for_effect(
                    n,
                    sd,
                    maf,
                    effect,
                )

                out.write(
                    f"{n}\t{maf:.2f}\t{effect:.1f}\t"
                    f"{sd:.6f}\t{ALPHA:.2e}\t"
                    f"{power:.8g}\n"
                )

def main():
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    verify_files()

    genotype_n = count_lines(INPUT_PREFIX.with_suffix(".fam"))

    create_keep_file()

    keep_n = count_lines(KEEP_FILE)

    phenotype_n, mean, sd = calculate_phenotype_stats()

    create_genotype_subset()

    subset_n = count_lines(GENO_PREFIX.with_suffix(".fam"))

    calculate_maf()

    variant_n, min_maf, median_maf, max_maf = get_maf_summary()

    write_power_table(phenotype_n, sd)

    print()
    print("Verification")
    print(f"Genotype samples: {genotype_n}")
    print(f"Phenotype samples: {phenotype_n}")
    print(f"Keep-file samples: {keep_n}")
    print(f"Subset samples: {subset_n}")
    print(f"Mean lifespan: {mean:.4f}")
    print(f"SD lifespan: {sd:.6f}")
    print(f"Variants: {variant_n}")
    print(f"Min MAF: {min_maf:.6f}")
    print(f"Median MAF: {median_maf:.6f}")
    print(f"Max MAF: {max_maf:.6f}")
    print(f"Power results: {POWER_FILE}")
    print()

if __name__ == "__main__":
    main()