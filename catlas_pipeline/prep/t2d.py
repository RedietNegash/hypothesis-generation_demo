import pandas as pd

GWAS_URL = "https://diagram-consortium.org/downloads.html"  # TODO: replace with the
# specific European-ancestry
# DIAMANTE file link once chosen

# Best-guess column map -> canonical schema. VERIFY against the real header
# (prepare_t2d prints the actual columns and raises if this map doesn't match).
COLUMN_MAP = {
    "chr": "Chr",
    "pos": "Pos",
    "ref": "EA",  # effect allele
    "alt": "NEA",  # non-effect allele
    "beta": "Beta",
    "pval": "Pvalue",
    "rsid": "SNP",
}

PVAL_THRESHOLD = 5e-8


def prepare_t2d(gwas_path: str, pval_threshold: float = PVAL_THRESHOLD) -> pd.DataFrame:
    """Returns canonical schema: variant | chr | pos | ref | alt | beta | pval | rsid"""
    header = pd.read_csv(gwas_path, sep="\t", compression="infer", nrows=0)
    missing = [c for c in COLUMN_MAP.values() if c not in header.columns]
    if missing:
        raise ValueError(
            f"[T2D] COLUMN_MAP expects {missing} but the real file has columns: "
            f"{header.columns.tolist()}. Update COLUMN_MAP in prep/t2d.py to match, "
            f"then re-run -- do not guess a mapping and proceed."
        )

    df = pd.read_csv(
        gwas_path,
        sep="\t",
        compression="infer",
        usecols=list(COLUMN_MAP.values()),
        dtype={COLUMN_MAP["chr"]: str},
    )
    df = df.rename(columns={v: k for k, v in COLUMN_MAP.items()})
    df = df[df["pval"] < pval_threshold].copy()

    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df["variant"] = (
        df["chr"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    )
    df = df[["variant", "chr", "pos", "ref", "alt", "beta", "pval", "rsid"]]
    print(f"[T2D] {len(df):,} variants below pval<{pval_threshold}")
    return df
