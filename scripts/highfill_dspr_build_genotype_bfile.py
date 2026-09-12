import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

BASE_DIR = Path(__file__).resolve().parent.parent

RAW_GENOTYPE_FILE = BASE_DIR / "data" / "Highfill_803RILs_genotype.txt"
LIFTOVER_CHAIN = BASE_DIR / "data" / "dm3ToDm6.over.chain.gz"

OUT_DIR = BASE_DIR / "data" / "highfill_finemap"
PED_FILE = OUT_DIR / "highfill_genotypes.ped"
MAP_FILE = OUT_DIR / "highfill_genotypes.map"
BFILE_RAW = OUT_DIR / "highfill_genotypes_raw"
A1_FILE = OUT_DIR / "highfill_a1_alleles.txt"
BFILE = OUT_DIR / "highfill_genotypes"

CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "23"}

HOM_A1_THRESHOLD = 0.8
HOM_A2_THRESHOLD = 0.2


def load_raw_genotypes() -> pd.DataFrame:
    df = pd.read_csv(
        RAW_GENOTYPE_FILE, sep="\t", header=None,
        names=["CHR", "POS", "RIL", "A1", "A2", "count1", "count2"],
    )
    df["SNP"] = df["CHR"] + ":" + df["POS"].astype(str)
    total = df["count1"] + df["count2"]
    frac_a1 = df["count1"] / total.replace(0, np.nan)

    df["call"] = np.select(
        [frac_a1 >= HOM_A1_THRESHOLD, frac_a1 <= HOM_A2_THRESHOLD],
        [2, 0],
        default=np.nan,
    )
    return df


def build_snp_table(df: pd.DataFrame) -> pd.DataFrame:
    snps = df.groupby("SNP", sort=False).first()[["CHR", "POS", "A1", "A2"]].reset_index()

    lo = LiftOver(str(LIFTOVER_CHAIN))
    lifted = []
    for chrom, pos in zip(snps["CHR"], snps["POS"]):
        r = lo.convert_coordinate(f"chr{chrom}", int(pos))
        lifted.append(r[0][1] if r else None)
    snps["POS_dm6"] = lifted
    snps = snps.dropna(subset=["POS_dm6"]).copy()
    snps["POS_dm6"] = snps["POS_dm6"].astype(int)
    snps["CHR_NUM"] = snps["CHR"].map(CHROM_MAP)

    print(f"  {len(snps):,} SNPs retained after dm3->dm6 liftover")
    return snps


def write_ped_map(df: pd.DataFrame, snps: pd.DataFrame) -> pd.DataFrame:
    snps = snps.copy()
    snps["_chr_sort"] = snps["CHR_NUM"].astype(int)
    snps = snps.sort_values(["_chr_sort", "POS_dm6"])
    snp_order = snps["SNP"].tolist()

    if PED_FILE.exists() and MAP_FILE.exists():
        print(f"PED/MAP already exist: {PED_FILE}, {MAP_FILE}")
        return snps

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(MAP_FILE, "w") as f:
        for _, row in snps.iterrows():
            f.write(f"{row['CHR_NUM']}\t{row['SNP']}\t0\t{row['POS_dm6']}\n")
    print(f"  MAP written: {MAP_FILE} ({len(snps):,} SNPs)")

    hom_a1 = (snps["A1"] + " " + snps["A1"]).to_numpy()
    hom_a2 = (snps["A2"] + " " + snps["A2"]).to_numpy()
    missing = np.full(len(snps), "0 0")

    calls = df.pivot(index="RIL", columns="SNP", values="call").reindex(columns=snp_order)
    call_matrix = calls.to_numpy()

    n_rils = 0
    with open(PED_FILE, "w") as f:
        for ril, call_row in zip(calls.index, call_matrix):
            genos = np.select([call_row == 2, call_row == 0], [hom_a1, hom_a2], default=missing)
            f.write(f"0\t{ril}\t0\t0\t0\t-9\t" + "\t".join(genos) + "\n")
            n_rils += 1
    print(f"  PED written: {PED_FILE} ({n_rils:,} RILs)")
    return snps


def make_bed(snps: pd.DataFrame):
    if BFILE.with_suffix(".bed").exists():
        print(f"bfile already exists: {BFILE}")
        return

    if not BFILE_RAW.with_suffix(".bed").exists():
        subprocess.run([
            "plink",
            "--file", str(PED_FILE.with_suffix("")),
            "--make-bed",
            "--out", str(BFILE_RAW),
        ], check=True)

    # PLINK's PED/MAP autoconversion has no real notion of which allele is A1
    # (unlike a VCF's REF/ALT) -- it just assigns A1 by scan order, which does
    # NOT reliably match the A1 we intended (the allele genotype.txt's dosage
    # counts, which our GWAS effect sizes are computed against). Force the
    # correct A1 explicitly via --a1-allele so downstream effect directions
    # (COJO, etc.) line up with the GWAS summary stats.
    snps[["SNP", "A1"]].to_csv(A1_FILE, sep="\t", header=False, index=False)
    subprocess.run([
        "plink",
        "--bfile", str(BFILE_RAW),
        "--a1-allele", str(A1_FILE),
        "--make-bed",
        "--out", str(BFILE),
    ], check=True)


if __name__ == "__main__":
    df = load_raw_genotypes()
    snps = build_snp_table(df)
    snps = write_ped_map(df, snps)
    make_bed(snps)
