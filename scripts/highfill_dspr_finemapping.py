import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy.stats import t as tdist

BASE_DIR = Path(__file__).resolve().parent.parent

GWAS_FILE = BASE_DIR / "data" / "Highfill_803RILs_GWAS_PCA.txt"
LIFTOVER_CHAIN = BASE_DIR / "data" / "dm3ToDm6.over.chain.gz"
LD_REF_BFILE = BASE_DIR / "data" / "highfill_finemap" / "highfill_genotypes"

OUT_DIR = BASE_DIR / "data" / "highfill_finemap"
FREQ_PREFIX = OUT_DIR / "highfill_genotypes_freq"
SIG_SNP_FILE = OUT_DIR / "highfill_significant_snps.tsv"
COJO_INPUT_FILE = OUT_DIR / "highfill_cojo_input.txt"
COJO_OUT_PREFIX = OUT_DIR / "cojo" / "highfill_lifespan_cojo"
REGIONS_DIR = OUT_DIR / "regions"

WINDOW_BP = 500_000  # DSPR RIL panels have far longer-range LD than DGRP (few
                      # recombination breakpoints per line), so a wider window
                      # than the 100kb used for the DGRP fly panel is warranted.

SIG_P_THRESHOLD = 1e-5  # genome-wide (5e-8) leaves 0 SNPs at this sample size
MIN_N = 100

GCTA_BIN = shutil.which("gcta64") or "/home/icog-bioai2/bin/gcta64"
PLINK_BIN = shutil.which("plink") or "/usr/local/bin/plink"


def build_snp_table() -> pd.DataFrame:
    """Merge the GWAS results with dm6-lifted positions and derive SE from BETA/P/N."""
    df = pd.read_csv(GWAS_FILE, sep="\t")

    lo = LiftOver(str(LIFTOVER_CHAIN))
    lifted = []
    for chrom, pos in zip(df["CHR"], df["POS"]):
        r = lo.convert_coordinate(f"chr{chrom}", int(pos))
        lifted.append(r[0][1] if r else None)
    df["POS_dm6"] = lifted
    df = df.dropna(subset=["POS_dm6"]).copy()
    df["POS_dm6"] = df["POS_dm6"].astype(int)

    df_p = df["P"].clip(lower=1e-300)
    t_stat = np.sign(df["BETA"]) * tdist.ppf(1 - df_p / 2, df=df["N"] - 2)
    df["se"] = df["BETA"] / t_stat

    print(f"Merged GWAS: {len(df):,} SNPs with dm6 positions and derived SE")
    return df


def compute_freq() -> pd.DataFrame:
    frq_file = FREQ_PREFIX.with_suffix(".frq")
    if not frq_file.exists():
        subprocess.run([
            PLINK_BIN,
            "--bfile", str(LD_REF_BFILE),
            "--freq",
            "--out", str(FREQ_PREFIX),
        ], check=True)

    freq = pd.read_csv(frq_file, sep=r"\s+")
    return freq[["SNP", "A1", "A2", "MAF"]]


def filter_significant_snps(gwas: pd.DataFrame, freq: pd.DataFrame) -> pd.DataFrame:
    if SIG_SNP_FILE.exists():
        print(f"Significant SNPs already extracted: {SIG_SNP_FILE}")
        return pd.read_csv(SIG_SNP_FILE, sep="\t")

    merged = gwas.merge(freq, on="SNP", how="inner")
    sig = merged[
        (merged["N"] >= MIN_N) & (merged["P"] <= SIG_P_THRESHOLD)
    ].sort_values("P").copy()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sig.to_csv(SIG_SNP_FILE, sep="\t", index=False)

    print(f"SNPs at p <= 5e-8 (standard genome-wide) : {int((gwas['P'] <= 5e-8).sum())}")
    print(f"SNPs at p <= {SIG_P_THRESHOLD:g} (suggestive, used here) : {len(sig)}")
    print(f"Saved: {SIG_SNP_FILE}")
    return sig


def write_cojo_input(sig: pd.DataFrame) -> None:
    cojo = sig.rename(columns={"MAF": "freq", "BETA": "b"})
    cojo = cojo[["SNP", "A1", "A2", "freq", "b", "se", "P", "N"]].rename(columns={"P": "p"})
    cojo["N"] = cojo["N"].astype(int)
    cojo.to_csv(COJO_INPUT_FILE, sep=" ", index=False)
    print(f"COJO input written: {COJO_INPUT_FILE} ({len(cojo)} SNPs)")


def run_cojo() -> None:
    jma_file = COJO_OUT_PREFIX.with_suffix(".jma.cojo")
    if jma_file.exists():
        print(f"COJO already run: {jma_file}")
        return

    COJO_OUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        GCTA_BIN,
        "--bfile", str(LD_REF_BFILE),
        "--cojo-file", str(COJO_INPUT_FILE),
        "--cojo-slct",
        "--cojo-p", str(SIG_P_THRESHOLD),
        # DSPR RILs have long-range LD (few recombination breakpoints per
        # line): the default 0.9 collinearity cutoff let 3 SNPs with pairwise
        # r=0.78-0.84 all into the joint model, producing wildly inflated,
        # unstable joint effect sizes (b~50 -> bJ~300). A stricter cutoff
        # correctly collapses them to the single strongest signal.
        "--cojo-collinear", "0.5",
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
    chrom_from_snp = signals["SNP"].str.split(":").str[0]

    for chrom, pos, snp in zip(chrom_from_snp, signals["bp"], signals["SNP"]):
        region = gwas[
            (gwas["CHR"] == chrom)
            & (gwas["POS_dm6"] >= pos - WINDOW_BP)
            & (gwas["POS_dm6"] <= pos + WINDOW_BP)
        ]
        out_file = REGIONS_DIR / f"chr{chrom}_pos{pos}_snps.tsv"
        region.to_csv(out_file, sep="\t", index=False)
        print(f"  {snp}: {len(region)} SNPs in +/-{WINDOW_BP // 1000}kb window -> {out_file}")


if __name__ == "__main__":
    gwas = build_snp_table()
    freq = compute_freq()
    sig = filter_significant_snps(gwas, freq)
    write_cojo_input(sig)
    run_cojo()
    signals = load_cojo_signals()
    extract_regions(gwas, signals)
