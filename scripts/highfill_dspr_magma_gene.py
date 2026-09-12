import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

BASE_DIR = Path(__file__).resolve().parent.parent
MAGMA_BIN = BASE_DIR / "tools" / "magma" / "magma"

BFILE = BASE_DIR / "data" / "Highfill_803RILs_MAGMA"
GWAS_PCA_FILE = BASE_DIR / "data" / "Highfill_803RILs_GWAS_PCA.txt"
GENE_LOC_FILE = BASE_DIR / "data" / "magma" / "gene_loc.txt"

MAGMA_DIR = BASE_DIR / "data" / "magma"
SNP_LOC_FILE = MAGMA_DIR / "highfill_snp_loc.txt"
PVAL_FILE = MAGMA_DIR / "highfill_pval.txt"

ANNOT_PREFIX = MAGMA_DIR / "highfill_lifespan"
ANNOT_FILE = MAGMA_DIR / "highfill_lifespan.genes.annot"

GENE_OUT_PREFIX = MAGMA_DIR / "highfill_lifespan_gene"
GENE_OUT_FILE = MAGMA_DIR / "highfill_lifespan_gene.genes.out"
GENE_RAW_FILE = MAGMA_DIR / "highfill_lifespan_gene.genes.raw"

GENE_SETS_FILE = MAGMA_DIR / "gene_sets.txt"
GENESET_OUT_PREFIX = MAGMA_DIR / "highfill_lifespan_geneset"
GENESET_OUT_FILE = MAGMA_DIR / "highfill_lifespan_geneset.gsa.out"

# Highfill/DSPR SNP IDs use fly chromosome arms (2L/2R/3L/3R/X); gene_loc.txt
# was built for the DGRP pipeline and numbers arms as plain integers instead
# (2L->1, 2R->2, 3L->3, 3R->4, 4->5, X stays X) -- same remap must be applied
# here or MAGMA can't match SNPs to genes by chromosome.
CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "X"}

# The Highfill/DSPR SNP positions (both the SNP ID text and the GWAS file's
# own POS column) are on the original dm3 assembly the DSPR panel was
# released on in 2012 -- gene_loc.txt is dm6-based, so positions must be
# lifted over before they can be matched to gene coordinates.
LIFTOVER_CHAIN = BASE_DIR / "data" / "dm3ToDm6.over.chain.gz"


def liftover_positions(df: pd.DataFrame) -> pd.DataFrame:
    lo = LiftOver(str(LIFTOVER_CHAIN))
    lifted_pos = []
    for chrom, pos in zip(df["CHR"], df["POS"]):
        result = lo.convert_coordinate(f"chr{chrom}", int(pos))
        lifted_pos.append(result[0][1] if result else None)
    df = df.copy()
    df["POS"] = lifted_pos
    n_dropped = df["POS"].isna().sum()
    if n_dropped:
        print(f"  {n_dropped:,} SNPs dropped (no dm3->dm6 liftover match)")
    df = df.dropna(subset=["POS"])
    df["POS"] = df["POS"].astype(int)
    return df


def build_snp_loc_file():
    if SNP_LOC_FILE.exists():
        print(f"SNP location file already exists: {SNP_LOC_FILE}")
        return

    df = pd.read_csv(GWAS_PCA_FILE, sep="\t", usecols=["SNP", "CHR", "POS"])
    df = liftover_positions(df)
    df["CHR"] = df["CHR"].map(CHROM_MAP)
    df[["SNP", "CHR", "POS"]].to_csv(SNP_LOC_FILE, sep="\t", header=False, index=False)
    print(f"  {len(df):,} SNP locations written to {SNP_LOC_FILE}")


def build_pval_file():
    if PVAL_FILE.exists():
        print(f"P-value file already exists: {PVAL_FILE}")
        return

    df = pd.read_csv(GWAS_PCA_FILE, sep="\t", usecols=["SNP", "P", "N"])
    df = df[["SNP", "P", "N"]]
    df.to_csv(PVAL_FILE, sep="\t", index=False)
    print(f"  {len(df):,} SNP p-values written to {PVAL_FILE}")


def run_annotate():
    if ANNOT_FILE.exists():
        print(f"Annotation file already exists: {ANNOT_FILE}")
        return

    subprocess.run([
        str(MAGMA_BIN),
        "--annotate", "nonhuman",
        "--snp-loc", str(SNP_LOC_FILE),
        "--gene-loc", str(GENE_LOC_FILE),
        "--out", str(ANNOT_PREFIX),
    ], check=True)


def run_gene_analysis():
    if GENE_OUT_FILE.exists():
        print(f"Gene analysis already exists: {GENE_OUT_FILE}")
        return

    subprocess.run([
        str(MAGMA_BIN),
        "--bfile", str(BFILE),
        "--pval", str(PVAL_FILE), "use=1,2", "ncol=N",
        "--gene-annot", str(ANNOT_FILE),
        "--out", str(GENE_OUT_PREFIX),
    ], check=True)


def run_geneset_analysis():
    if GENESET_OUT_FILE.exists():
        print(f"Gene-set results already exist: {GENESET_OUT_FILE}")
        return

    subprocess.run([
        str(MAGMA_BIN),
        "--gene-results", str(GENE_RAW_FILE),
        "--set-annot", str(GENE_SETS_FILE),
        "--out", str(GENESET_OUT_PREFIX),
    ], check=True)


def report_geneset_results():
    df = pd.read_csv(GENESET_OUT_FILE, sep=r"\s+", comment="#")
    n = len(df)
    bonf = 0.05 / n

    df = df.sort_values("P").reset_index(drop=True)
    rank = np.arange(1, n + 1)
    df["BH_q"] = df["P"] * n / rank
    df["BH_q"] = df["BH_q"][::-1].cummin()[::-1]

    print(f"\nTotal gene sets tested: {n}")
    print(f"Bonferroni (0.05/{n} = {bonf:.2e}): {(df['P'] < bonf).sum()} sets")
    print(f"Nominal P<0.05: {(df['P'] < 0.05).sum()} sets")
    print(f"BH-FDR<0.05: {(df['BH_q'] < 0.05).sum()} sets")
    print(f"BH-FDR<0.20: {(df['BH_q'] < 0.20).sum()} sets")

    print("\nTop 15 gene sets:")
    print(df[["FULL_NAME", "NGENES", "BETA", "SE", "P", "BH_q"]].head(15).to_string(index=False))


def report_results():
    df = pd.read_csv(GENE_OUT_FILE, sep=r"\s+", comment="#")
    n = len(df)
    bonf = 0.05 / n

    df = df.sort_values("P").reset_index(drop=True)
    rank = np.arange(1, n + 1)
    df["BH_q"] = df["P"] * n / rank
    df["BH_q"] = df["BH_q"][::-1].cummin()[::-1]

    print(f"\nTotal genes tested: {n}")
    print(f"Bonferroni (0.05/{n} = {bonf:.2e}): {(df['P'] < bonf).sum()} genes")
    print(f"Nominal P<0.05: {(df['P'] < 0.05).sum()} genes")
    print(f"BH-FDR<0.05: {(df['BH_q'] < 0.05).sum()} genes")
    print(f"BH-FDR<0.10: {(df['BH_q'] < 0.10).sum()} genes")

    print("\nTop 20 genes:")
    print(df[["GENE", "CHR", "START", "STOP", "NSNPS", "N", "ZSTAT", "P", "BH_q"]].head(20).to_string(index=False))


if __name__ == "__main__":
    build_snp_loc_file()
    build_pval_file()
    run_annotate()
    run_gene_analysis()
    report_results()
    run_geneset_analysis()
    report_geneset_results()
