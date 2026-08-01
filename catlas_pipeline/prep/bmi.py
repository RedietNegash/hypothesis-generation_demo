import pandas as pd

GWAS_URL = (
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/"
    "additive-tsvs/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz"
)
VAR_URL = (
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/variants.tsv.bgz"
)


GWAS_COLUMN_MAP = {
    "variant": "variant",
    "beta": "beta",
    "pval": "pval",
}
VAR_COLUMN_MAP = {
    "variant": "variant",
    "rsid": "rsid",
    "chr": "chr",
    "pos": "pos",
    "ref": "ref",
    "alt": "alt",
}

PVAL_THRESHOLD = 5e-8


def prepare_bmi(
    gwas_path: str, var_path: str, pval_threshold: float = PVAL_THRESHOLD
) -> pd.DataFrame:
    """Returns canonical schema: variant | chr | pos | ref | alt | beta | pval | rsid"""
    gwas = pd.read_csv(
        gwas_path, sep="\t", compression="gzip", usecols=list(GWAS_COLUMN_MAP.values())
    )
    gwas = gwas.rename(columns={v: k for k, v in GWAS_COLUMN_MAP.items()})
    gwas = gwas[gwas["pval"] < pval_threshold].copy()
    print(f"[BMI] {len(gwas):,} variants below pval<{pval_threshold}")

    var = pd.read_csv(
        var_path,
        sep="\t",
        compression="gzip",
        usecols=list(VAR_COLUMN_MAP.values()),
        dtype={"chr": str},
    )
    var = var.rename(columns={v: k for k, v in VAR_COLUMN_MAP.items()})

    df = gwas.merge(var, on="variant", how="inner")
    df = df[["variant", "chr", "pos", "ref", "alt", "beta", "pval", "rsid"]]
    print(f"[BMI] {len(df):,} variants after joining rsID annotation")
    return df
