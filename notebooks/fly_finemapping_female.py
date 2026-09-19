#!/usr/bin/env python3
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
FLY_CHROMS = ("2L", "2R", "3L", "3R", "4", "X")

PHENO_NAME = "S18_1537_F"
GLM_DIR = BASE_DIR / "data" / "gwas" / "tmp"

MERGED_QC_SOURCE = GLM_DIR / "merged_qc"

OUT_DIR = BASE_DIR / "data" / "finemap" / "female"
SIG_SNP_FILE = OUT_DIR / "female_significant_snps.tsv"
COJO_INPUT_FILE = OUT_DIR / "female_cojo_input.txt"
COJO_OUT_PREFIX = OUT_DIR / "cojo" / "female_lifespan_cojo"
REGIONS_DIR = OUT_DIR / "regions"

CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "23", "23": "23"}
LD_REF_BFILE = OUT_DIR / "bfile" / "merged_qc_numeric"

WINDOW_BP = 100_000

SIG_P_THRESHOLD = 1e-5
MIN_MAF = 0.05
MIN_N = 100

GCTA_BIN = "gcta64"

GWAS_COLUMN_MAP = {
    "#CHROM": "CHR",
    "ID": "SNP",
    "A1_FREQ": "freq",
    "BETA": "b",
    "SE": "se",
    "P": "p",
    "OBS_CT": "N",
    "OMITTED": "A2",
}
GWAS_COLUMNS = ["CHR", "POS", "SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
NUMERIC_COLUMNS = ["POS", "freq", "b", "se", "p", "N"]
COJO_COLUMNS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
COJO_RESULT_COLUMNS = ["Chr", "SNP", "bp", "b", "p", "bJ", "pJ"]
COJO_RESULT_NUMERIC_COLUMNS = ["bp", "b", "p", "bJ", "pJ"]


def merge_gwas_sumstats() -> pd.DataFrame:
    frames = []
    for chrom in FLY_CHROMS:
        input_file = GLM_DIR / f"lifespan_{chrom}.{PHENO_NAME}.glm.linear"
        if not input_file.is_file():
            raise FileNotFoundError(f"Missing female GWAS output for {chrom}: {input_file}")

        chromosome_gwas = pd.read_csv(input_file, sep="\t", low_memory=False)
        required_columns = set(GWAS_COLUMN_MAP) | {"POS", "A1"}
        missing_columns = sorted(required_columns - set(chromosome_gwas.columns))
        if missing_columns:
            raise ValueError(f"{input_file} is missing columns: {', '.join(missing_columns)}")

        chromosome_gwas = chromosome_gwas.rename(columns=GWAS_COLUMN_MAP)
        frames.append(chromosome_gwas[GWAS_COLUMNS])

    gwas = pd.concat(frames, ignore_index=True)
    for column in NUMERIC_COLUMNS:
        gwas[column] = pd.to_numeric(gwas[column], errors="coerce")

    row_count = len(gwas)
    gwas = gwas.dropna(subset=NUMERIC_COLUMNS).copy()
    gwas["CHR"] = gwas["CHR"].astype(str)
    dropped_count = row_count - len(gwas)

    print(f"Merged female GWAS: {len(gwas):,} SNPs across {len(FLY_CHROMS)} chromosomes")
    if dropped_count:
        print(f"Dropped {dropped_count:,} rows with missing or invalid numeric values")
    return gwas


def filter_significant_snps(gwas: pd.DataFrame) -> pd.DataFrame:
    maf = np.minimum(gwas["freq"], 1 - gwas["freq"])
    valid = (
        gwas[["SNP", "A1", "A2"]].notna().all(axis=1)
        & gwas["freq"].between(0, 1)
        & gwas["p"].between(0, 1)
        & (gwas["se"] > 0)
        & (gwas["N"] >= MIN_N)
    )
    eligible = gwas.loc[valid & (maf >= MIN_MAF)].copy()
    significant = eligible.loc[eligible["p"] <= SIG_P_THRESHOLD].sort_values("p")

    if significant.empty:
        raise ValueError(
            f"No SNPs passed MAF >= {MIN_MAF:g}, N >= {MIN_N}, "
            f"and P <= {SIG_P_THRESHOLD:g}"
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    significant.to_csv(SIG_SNP_FILE, sep="\t", index=False)

    genome_wide_count = int((eligible["p"] <= 5e-8).sum())
    print(f"Eligible SNPs after MAF and sample-size filters: {len(eligible):,}")
    print(f"SNPs at P <= 5e-8: {genome_wide_count:,}")
    print(f"SNPs at P <= {SIG_P_THRESHOLD:g}: {len(significant):,}")
    print(f"Saved significant SNPs: {SIG_SNP_FILE}")
    return significant


def write_cojo_input(significant_snps: pd.DataFrame) -> None:
    duplicate_snps = significant_snps.loc[
        significant_snps["SNP"].duplicated(keep=False), "SNP"
    ].unique()
    if len(duplicate_snps):
        duplicate_list = ", ".join(map(str, duplicate_snps[:5]))
        raise ValueError(f"COJO input contains duplicate SNP IDs: {duplicate_list}")

    cojo = significant_snps[COJO_COLUMNS].copy()
    rounded_sample_size = np.rint(cojo["N"])
    if not np.allclose(cojo["N"], rounded_sample_size):
        raise ValueError("COJO sample sizes must be whole numbers")

    cojo["N"] = rounded_sample_size.astype(int)
    COJO_INPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    cojo.to_csv(COJO_INPUT_FILE, sep=" ", index=False)
    print(f"Saved COJO input: {COJO_INPUT_FILE} ({len(cojo):,} SNPs)")


def prepare_cojo_bfile() -> None:
    source_files = [MERGED_QC_SOURCE.with_suffix(ext) for ext in (".bed", ".bim", ".fam")]
    missing_files = [path for path in source_files if not path.is_file()]
    if missing_files:
        missing_list = "\n".join(f"  {path}" for path in missing_files)
        raise FileNotFoundError(f"Missing PLINK LD-reference files:\n{missing_list}")

    LD_REF_BFILE.parent.mkdir(parents=True, exist_ok=True)
    for extension in (".bed", ".fam"):
        source = MERGED_QC_SOURCE.with_suffix(extension).resolve()
        destination = LD_REF_BFILE.with_suffix(extension)
        if destination.is_symlink() and destination.resolve() != source:
            destination.unlink()
        elif destination.exists() and not destination.is_symlink():
            raise FileExistsError(f"Refusing to replace existing file: {destination}")
        if not destination.exists():
            destination.symlink_to(source)

    source_bim = MERGED_QC_SOURCE.with_suffix(".bim")
    target_bim = LD_REF_BFILE.with_suffix(".bim")
    temporary_bim = target_bim.with_suffix(".bim.tmp")
    variant_count = 0

    try:
        with source_bim.open(encoding="utf-8") as input_file, temporary_bim.open(
            "w", encoding="utf-8"
        ) as output_file:
            for line_number, line in enumerate(input_file, start=1):
                fields = line.split()
                if len(fields) != 6:
                    raise ValueError(
                        f"Malformed BIM row {line_number} in {source_bim}: "
                        f"expected 6 fields, found {len(fields)}"
                    )

                chrom, snp, cm, pos, allele_1, allele_2 = fields
                if chrom not in CHROM_MAP:
                    raise ValueError(f"Unsupported chromosome {chrom!r} in {source_bim}")

                output_file.write(
                    f"{CHROM_MAP[chrom]}\t{snp}\t{cm}\t{pos}\t{allele_1}\t{allele_2}\n"
                )
                variant_count += 1

        temporary_bim.replace(target_bim)
    finally:
        temporary_bim.unlink(missing_ok=True)

    print(f"Prepared COJO LD reference: {LD_REF_BFILE} ({variant_count:,} SNPs)")


def run_cojo() -> None:
    jma_file = COJO_OUT_PREFIX.with_suffix(".jma.cojo")
    if jma_file.is_file():
        print(f"Using existing COJO result: {jma_file}")
        return

    gcta_executable = shutil.which(GCTA_BIN)
    if gcta_executable is None:
        gcta_path = Path(GCTA_BIN).expanduser()
        if not gcta_path.is_file():
            raise FileNotFoundError(
                f"Could not find GCTA executable {GCTA_BIN!r}; add gcta64 to PATH"
            )
        gcta_executable = str(gcta_path.resolve())

    required_inputs = [
        COJO_INPUT_FILE,
        LD_REF_BFILE.with_suffix(".bed"),
        LD_REF_BFILE.with_suffix(".bim"),
        LD_REF_BFILE.with_suffix(".fam"),
    ]
    missing_inputs = [path for path in required_inputs if not path.is_file()]
    if missing_inputs:
        missing_list = "\n".join(f"  {path}" for path in missing_inputs)
        raise FileNotFoundError(f"Missing COJO input files:\n{missing_list}")

    COJO_OUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    command = [
        gcta_executable,
        "--bfile",
        str(LD_REF_BFILE),
        "--maf",
        str(MIN_MAF),
        "--cojo-file",
        str(COJO_INPUT_FILE),
        "--cojo-slct",
        "--cojo-p",
        str(SIG_P_THRESHOLD),
        "--out",
        str(COJO_OUT_PREFIX),
    ]
    print(f"Running: {shlex.join(command)}")
    subprocess.run(command, check=True)

    if not jma_file.is_file():
        raise RuntimeError(f"GCTA finished without creating the expected result: {jma_file}")
    print(f"Saved COJO result: {jma_file}")


def load_cojo_signals() -> pd.DataFrame:
    jma_file = COJO_OUT_PREFIX.with_suffix(".jma.cojo")
    if not jma_file.is_file():
        raise FileNotFoundError(f"Missing COJO result: {jma_file}")

    signals = pd.read_csv(jma_file, sep=r"\s+")
    missing_columns = sorted(set(COJO_RESULT_COLUMNS) - set(signals.columns))
    if missing_columns:
        raise ValueError(f"{jma_file} is missing columns: {', '.join(missing_columns)}")
    if signals.empty:
        raise ValueError(f"COJO selected no independent signals in {jma_file}")

    for column in COJO_RESULT_NUMERIC_COLUMNS:
        signals[column] = pd.to_numeric(signals[column], errors="coerce")
    if signals[COJO_RESULT_NUMERIC_COLUMNS].isna().any(axis=None):
        raise ValueError(f"COJO result contains missing or invalid numeric values: {jma_file}")
    if not signals["p"].between(0, 1).all() or not signals["pJ"].between(0, 1).all():
        raise ValueError(f"COJO result contains P-values outside [0, 1]: {jma_file}")

    rounded_positions = np.rint(signals["bp"])
    if not np.allclose(signals["bp"], rounded_positions) or (rounded_positions < 1).any():
        raise ValueError(f"COJO result contains invalid base-pair positions: {jma_file}")
    signals["bp"] = rounded_positions.astype(int)

    duplicate_snps = signals.loc[signals["SNP"].duplicated(keep=False), "SNP"].unique()
    if len(duplicate_snps):
        duplicate_list = ", ".join(map(str, duplicate_snps[:5]))
        raise ValueError(f"COJO result contains duplicate SNP IDs: {duplicate_list}")

    signals = signals.sort_values("pJ").reset_index(drop=True)
    print(f"\n{len(signals)} independent signal(s) selected by COJO:")
    print(signals[COJO_RESULT_COLUMNS].to_string(index=False))
    return signals


def extract_regions(gwas: pd.DataFrame, signals: pd.DataFrame) -> None:
    REGIONS_DIR.mkdir(parents=True, exist_ok=True)
    chrom_from_snp = signals["SNP"].str.split("_").str[0]

    for chrom, pos, snp in zip(chrom_from_snp, signals["bp"], signals["SNP"]):
        region = gwas[
            (gwas["CHR"] == chrom)
            & (gwas["POS"] >= pos - WINDOW_BP)
            & (gwas["POS"] <= pos + WINDOW_BP)
        ]
        out_file = REGIONS_DIR / f"chr{chrom}_pos{pos}_snps.tsv"
        region.to_csv(out_file, sep="\t", index=False)
        print(f"  {snp}: {len(region)} SNPs in +/-{WINDOW_BP // 1000}kb window -> {out_file}")


if __name__ == "__main__":
    gwas = merge_gwas_sumstats()
    sig = filter_significant_snps(gwas)
    write_cojo_input(sig)
    prepare_cojo_bfile()
    run_cojo()
    signals = load_cojo_signals()
    extract_regions(gwas, signals)
