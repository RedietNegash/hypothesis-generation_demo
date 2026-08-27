import pandas as pd

GWAS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST007001-GCST008000/GCST007320/AD_sumstats_Jansenetal_2019sept.txt.gz"
)  # original ctg.cncr.nl host is dead; this is EBI GWAS Catalog's mirror (GCST007320)

COLUMN_MAP = {
    "chr": "CHR",
    "pos": "BP",
    "ref": "A1",  # effect allele
    "alt": "A2",
    "beta": "BETA",
    "pval": "P",
    "rsid": "SNP",
}

PVAL_THRESHOLD = 5e-8


def prepare_ad(gwas_path: str, pval_threshold: float = PVAL_THRESHOLD) -> pd.DataFrame:
    """Returns canonical schema: variant | chr | pos | ref | alt | beta | pval | rsid"""
    header = pd.read_csv(gwas_path, sep="\t", compression="infer", nrows=0)
    missing = [c for c in COLUMN_MAP.values() if c not in header.columns]
    if missing:
        raise ValueError(
            f"[AD] COLUMN_MAP expects {missing} but the real file has columns: "
            f"{header.columns.tolist()}. Update COLUMN_MAP in prep/ad.py to match, "
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
    print(f"[AD] {len(df):,} variants below pval<{pval_threshold}")
    return df
