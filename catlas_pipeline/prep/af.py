import pandas as pd

GWAS_URL = (
    "http://csg.sph.umich.edu/willer/public/afib2018/"
    "nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz"
)

# Column map -> canonical schema, verified against the actual file header:
# MarkerName rs_dbSNP147 CHR POS_GRCh37 A1 A2 Freq_A2 Effect_A2 StdErr Pvalue
COLUMN_MAP = {
    "chr": "CHR",
    "pos": "POS_GRCh37",
    "ref": "A1",
    "alt": "A2",
    "beta": "Effect_A2",
    "pval": "Pvalue",
    "rsid": "rs_dbSNP147",
}

PVAL_THRESHOLD = 5e-8


def prepare_af(gwas_path: str, pval_threshold: float = PVAL_THRESHOLD) -> pd.DataFrame:
    """Returns canonical schema: variant | chr | pos | ref | alt | beta | pval | rsid"""
    df = pd.read_csv(
        gwas_path,
        sep="\t",
        compression="gzip",
        usecols=list(COLUMN_MAP.values()),
        dtype={COLUMN_MAP["chr"]: str},
    )
    df = df.rename(columns={v: k for k, v in COLUMN_MAP.items()})
    df = df[df["pval"] < pval_threshold].copy()

    df["variant"] = (
        df["chr"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    )
    df = df[["variant", "chr", "pos", "ref", "alt", "beta", "pval", "rsid"]]
    print(f"[AF] {len(df):,} variants below pval<{pval_threshold}")
    return df
