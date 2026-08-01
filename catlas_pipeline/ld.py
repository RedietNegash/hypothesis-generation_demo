import subprocess
import pandas as pd


def load_reference_bim(bim_prefix: str, min_rows: int = 1_000_000) -> pd.DataFrame:
    """Load a PLINK .bim (LD reference panel, hg19/GRCh37)."""
    bim = pd.read_csv(
        f"{bim_prefix}.bim",
        sep="\t",
        header=None,
        names=["chr", "rsid", "cm", "pos", "a1", "a2"],
    )
    assert len(bim) > min_rows, (
        f"Only read {len(bim):,} rows from {bim_prefix}.bim -- reference panel may "
        f"still be syncing or the prefix is wrong, re-check and re-run"
    )
    print(f"Reference .bim: {len(bim):,} variants (hg19/GRCh37)")
    return bim


load_eur_bim = load_reference_bim


def backfill_rsids_from_bim(
    df: pd.DataFrame, bim: pd.DataFrame, label: str = ""
) -> pd.DataFrame:
    """Replace the ``rsid`` column by a chr+pos join against a reference .bim.

    Use when a GWAS file has no real dbSNP rsIDs (e.g. its identifier column is
    ``chr:pos``). df and bim must be the SAME genome build (the reference panel is
    hg19), so run this before liftover. Variants with no chr:pos match are dropped.
    """
    df = df.copy()
    bim = bim.copy()
    df["chr"] = df["chr"].astype(str)
    bim["chr"] = bim["chr"].astype(str)
    df["pos"] = df["pos"].astype("int64")
    bim["pos"] = bim["pos"].astype("int64")
    before = len(df)
    df = df.drop(columns=["rsid"]).merge(
        bim[["chr", "pos", "rsid"]].drop_duplicates(subset=["chr", "pos"]),
        on=["chr", "pos"],
        how="inner",
    )
    print(
        f"[{label}] rsID backfill via chr:pos join: {len(df):,}/{before:,} "
        f"variants matched to reference-panel rsIDs"
    )
    canonical = ["variant", "chr", "pos", "ref", "alt", "beta", "pval", "rsid"]
    return df[[c for c in canonical if c in df.columns]]


def match_rsids(df: pd.DataFrame, bim_rsids: set, label: str = "") -> pd.DataFrame:
    matched = df["rsid"].isin(bim_rsids)
    rate = matched.mean() * 100
    print(f"[{label}] rsID match rate: {rate:.1f}% ({matched.sum():,}/{len(df):,})")
    if rate < 5 and len(df):
        print(
            f"[{label}] WARNING: almost no rsIDs matched the reference panel. If this "
            f"GWAS lacks real dbSNP rsIDs (e.g. its ID column is chr:pos), re-run with "
            f"--backfill-rsids (CLI) or ld.backfill_rsids_from_bim(...) before matching."
        )
    return df[matched].copy()


def ld_prune(
    matched_df: pd.DataFrame,
    eur_bim_prefix: str,
    plink_bin: str,
    dataset_dir: str,
    prune_params: tuple,
    label: str = "",
) -> pd.DataFrame:

    rsid_file = f"{dataset_dir}/rsids_for_plink.txt"
    with open(rsid_file, "w") as f:
        for rsid in matched_df["rsid"]:
            f.write(f"{rsid}\n")

    out_prefix = f"{dataset_dir}/ld_pruned"
    cmd = [
        plink_bin,
        "--bfile",
        eur_bim_prefix,
        "--extract",
        rsid_file,
        "--indep-pairwise",
        *prune_params,
        "--out",
        out_prefix,
        "--memory",
        "8000",
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    with open(f"{out_prefix}.prune.in") as f:
        pruned_rsids = set(f.read().splitlines())

    df_pruned = matched_df[matched_df["rsid"].isin(pruned_rsids)].copy()
    print(
        f"[{label}] LD pruning ({'/'.join(prune_params)}): "
        f"{len(matched_df):,} -> {len(df_pruned):,} independent variants"
    )
    return df_pruned
