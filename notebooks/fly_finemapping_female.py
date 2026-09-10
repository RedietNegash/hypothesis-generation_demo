#!/usr/bin/env python3
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
FLY_CHROMS = ["2L", "2R", "3L", "3R", "4", "X"]

PHENO_NAME = "S18_1537_F"
GLM_DIR = BASE_DIR / "data" / "gwas" / "tmp"

MERGED_QC_SOURCE = GLM_DIR / "merged_qc"

OUT_DIR = BASE_DIR / "data" / "finemap" / "female"
SIG_SNP_FILE = OUT_DIR / "female_significant_snps.tsv"
COJO_INPUT_FILE = OUT_DIR / "female_cojo_input.txt"
COJO_OUT_PREFIX = OUT_DIR / "cojo" / "female_lifespan_cojo"
REGIONS_DIR = OUT_DIR / "regions"

CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "23": "23"}
LD_REF_BFILE = OUT_DIR / "bfile" / "merged_qc_numeric"

WINDOW_BP = 100_000

SIG_P_THRESHOLD = 1e-5
MIN_MAF = 0.05
MIN_N = 100

GCTA_BIN = shutil.which("gcta64") or "/home/icog-bioai2/bin/gcta64"


def merge_gwas_sumstats() -> pd.DataFrame:
    frames = []
    for chrom in FLY_CHROMS:
        f = GLM_DIR / f"lifespan_{chrom}.{PHENO_NAME}.glm.linear"
        if not f.exists():
            raise FileNotFoundError(f"Missing female GWAS output for {chrom}: {f}")
        frames.append(pd.read_csv(f, sep="\t"))
    df = pd.concat(frames, ignore_index=True)

    df = df.rename(columns={
        "#CHROM": "CHR", "ID": "SNP", "A1_FREQ": "freq",
        "BETA": "b", "SE": "se", "P": "p", "OBS_CT": "N", "OMITTED": "A2",
    })
    df = df[["CHR", "POS", "SNP", "A1", "A2", "freq", "b", "se", "p", "N"]].copy()

    for col in ["freq", "b", "se", "p", "N"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["freq", "b", "se", "p", "N"])

    print(f"Merged female GWAS: {len(df):,} SNPs across {len(FLY_CHROMS)} chromosomes")
    return df


def filter_significant_snps(gwas: pd.DataFrame) -> pd.DataFrame:
    if SIG_SNP_FILE.exists():
        print(f"Significant SNPs already extracted: {SIG_SNP_FILE}")
        return pd.read_csv(SIG_SNP_FILE, sep="\t")

    maf = np.minimum(gwas["freq"], 1 - gwas["freq"])
    sig = gwas[
        (maf > MIN_MAF) & (gwas["N"] >= MIN_N) & (gwas["p"] <= SIG_P_THRESHOLD)
    ].sort_values("p").copy()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sig.to_csv(SIG_SNP_FILE, sep="\t", index=False)

    n_gw = int((gwas["p"] <= 5e-8).sum())
    print(f"SNPs at p <= 5e-8 (standard genome-wide)      : {n_gw}")
    print(f"SNPs at p <= {SIG_P_THRESHOLD:g} (suggestive, used here) : {len(sig)}")
    print(f"Saved: {SIG_SNP_FILE}")
    return sig


def write_cojo_input(sig: pd.DataFrame) -> None:
    cojo = sig[["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]].copy()
    cojo["N"] = cojo["N"].astype(int)
    cojo.to_csv(COJO_INPUT_FILE, sep=" ", index=False)
    print(f"COJO input written: {COJO_INPUT_FILE} ({len(cojo)} SNPs)")


def prepare_cojo_bfile() -> None:
    bim_out = LD_REF_BFILE.with_suffix(".bim")
    if bim_out.exists():
        print(f"COJO bfile already prepared: {LD_REF_BFILE}")
        return

    LD_REF_BFILE.parent.mkdir(parents=True, exist_ok=True)
    for ext in [".bed", ".fam"]:
        dst = LD_REF_BFILE.with_suffix(ext)
        if not dst.exists():
            dst.symlink_to(MERGED_QC_SOURCE.with_suffix(ext).resolve())

    n = 0
    with open(MERGED_QC_SOURCE.with_suffix(".bim")) as fh, open(bim_out, "w") as out:
        for line in fh:
            chrom, snp, cm, pos, a1, a2 = line.rstrip("\n").split("\t")
            out.write(f"{CHROM_MAP[chrom]}\t{snp}\t{cm}\t{pos}\t{a1}\t{a2}\n")
            n += 1

    print(f"  {n:,} SNPs remapped into {bim_out}")


def run_cojo() -> None:
    jma_file = COJO_OUT_PREFIX.with_suffix(".jma.cojo")
    if jma_file.exists():
        print(f"COJO already run: {jma_file}")
        return

    COJO_OUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        GCTA_BIN,
        "--bfile", str(LD_REF_BFILE),
        "--maf", str(MIN_MAF),
        "--cojo-file", str(COJO_INPUT_FILE),
        "--cojo-slct",
        "--cojo-p", str(SIG_P_THRESHOLD),
        "--out", str(COJO_OUT_PREFIX),
    ], check=True)
    print(f"COJO finished: {jma_file}")


def load_cojo_signals() -> pd.DataFrame:
    jma_file = COJO_OUT_PREFIX.with_suffix(".jma.cojo")
    signals = pd.read_csv(jma_file, sep=r"\s+").sort_values("pJ").reset_index(drop=True)
    print(f"\n{len(signals)} independent signal(s) selected by COJO:")
    print(signals[["Chr", "SNP", "bp", "b", "p", "bJ", "pJ"]].to_string(index=False))
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
