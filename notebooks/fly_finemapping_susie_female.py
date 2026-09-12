#!/usr/bin/env python3
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")

REGIONS_DIR = BASE_DIR / "data" / "finemap" / "female" / "regions"
LD_REF_BFILE = BASE_DIR / "data" / "finemap" / "female" / "bfile" / "merged_qc_numeric"
WORK_DIR = BASE_DIR / "data" / "finemap" / "female" / "susie"
OUT_DIR = BASE_DIR / "data" / "finemap" / "female" / "susie_results"

PLINK_BIN = shutil.which("plink") or "/usr/local/bin/plink"

COVERAGE = 0.95
L = 10

susieR = importr("susieR")
numpy2ri.activate()
pandas2ri.activate()


def region_files():
    return sorted(REGIONS_DIR.glob("chr*_snps.tsv"))


def region_label(region_file: Path) -> str:
    return region_file.stem.replace("_snps", "")


def prepare_region_bfile(region_file: Path, label: str) -> Path:
    bfile_out = WORK_DIR / label
    if bfile_out.with_suffix(".bim").exists():
        return bfile_out

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    snplist_file = WORK_DIR / f"{label}.snplist"
    region = pd.read_csv(region_file, sep="\t")
    snplist_file.write_text("\n".join(region["SNP"]) + "\n")

    subprocess.run([
        PLINK_BIN,
        "--bfile", str(LD_REF_BFILE),
        "--extract", str(snplist_file),
        "--make-bed",
        "--out", str(bfile_out),
    ], check=True, capture_output=True, text=True)

    subprocess.run([
        PLINK_BIN,
        "--bfile", str(bfile_out),
        "--r", "square",
        "--out", str(bfile_out),
    ], check=True, capture_output=True, text=True)

    return bfile_out


def load_region_data(region_file: Path, bfile: Path):
    bim = pd.read_csv(
        bfile.with_suffix(".bim"), sep="\t", header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    ld = np.loadtxt(bfile.with_suffix(".ld"))

    region = pd.read_csv(region_file, sep="\t").set_index("SNP")
    ordered = bim.join(region[["b", "se", "p", "N"]], on="SNP")

    keep = ordered["b"].notna().to_numpy()
    ordered = ordered[keep].reset_index(drop=True)
    ld = ld[np.ix_(keep, keep)]

    keep_nan = ~np.isnan(ld).any(axis=1)
    ordered = ordered[keep_nan].reset_index(drop=True)
    ld = ld[np.ix_(keep_nan, keep_nan)]

    return ordered, ld


def run_susie_rss(ordered: pd.DataFrame, ld: np.ndarray):
    bhat = ro.FloatVector(ordered["b"].to_numpy())
    shat = ro.FloatVector(ordered["se"].to_numpy())
    n = int(ordered["N"].median())
    R = ro.r["matrix"](ro.FloatVector(ld.flatten(order="F")), nrow=ld.shape[0])

    fit = susieR.susie_rss(
        bhat=bhat, shat=shat, R=R, n=n, L=L,
        estimate_residual_variance=True, verbose=False,
    )
    pip = np.array(fit.rx2("pip"))

    cs_result = susieR.susie_get_cs(fit, coverage=COVERAGE, Xcorr=R)
    cs_list = cs_result.rx2("cs")

    cs_assignment = np.full(len(ordered), np.nan)
    if cs_list != ro.NULL:
        names = [str(n) for n in cs_list.names]
        for cs_idx, cs_name in enumerate(names, start=1):
            member_idx = np.array(cs_list.rx2(cs_name)) - 1
            cs_assignment[member_idx] = cs_idx

    return pip, cs_assignment


def finemap_region(region_file: Path) -> None:
    label = region_label(region_file)
    out_file = OUT_DIR / f"{label}_susie.tsv"
    if out_file.exists():
        print(f"Already fine-mapped: {out_file}")
        return

    bfile = prepare_region_bfile(region_file, label)
    ordered, ld = load_region_data(region_file, bfile)

    if len(ordered) < 2:
        print(f"{label}: too few matched SNPs ({len(ordered)}), skipping")
        return

    pip, cs = run_susie_rss(ordered, ld)
    ordered["PIP"] = pip
    ordered["CS"] = cs
    ordered = ordered.sort_values("PIP", ascending=False)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ordered.to_csv(out_file, sep="\t", index=False)

    n_cs = int(np.nanmax(cs)) if np.any(~np.isnan(cs)) else 0
    top = ordered.iloc[0]
    print(f"{label}: {len(ordered)} SNPs, {n_cs} credible set(s), top PIP {top['SNP']} = {top['PIP']:.3f}")
    print(f"  Saved: {out_file}")


if __name__ == "__main__":
    for region_file in region_files():
        finemap_region(region_file)
