import pandas as pd

CANONICAL = ["variant", "chr", "pos", "ref", "alt", "beta", "pval", "rsid"]


def prepare_gwas(
    gwas_path: str,
    column_map: dict,
    *,
    pval_threshold: float = 5e-8,
    sep: str = "\t",
    compression: str = "infer",
    uppercase_alleles: bool = False,
    label: str = "GWAS",
) -> pd.DataFrame:
    """Read an arbitrary GWAS summary-statistics file into the pipeline's canonical
    schema:  ``variant | chr | pos | ref | alt | beta | pval | rsid``.

    ``column_map`` maps a canonical field -> the source column name in the file.
      Required keys: ``chr``, ``pos``, ``ref``, ``alt``, ``pval``.
      Optional keys:
        ``beta`` -- effect size. If absent, filled with 0.0 (allele reorientation
                    then cannot flip the effect sign; fine if you only need which
                    cell types are enriched, not the direction of effect).
        ``rsid`` -- dbSNP id. If absent, a ``chr:pos`` placeholder is used; pair
                    with rsID backfill (``ld.backfill_rsids_from_bim``) so LD
                    pruning has real ids to work with.

    ``ref``/``alt`` are treated as the two alleles at the locus; which is the
    genome reference is sorted out later by ``genome.reorient_alleles``, so it is
    fine to map them to effect/non-effect allele columns in either order.
    """
    required = ["chr", "pos", "ref", "alt", "pval"]
    missing = [k for k in required if k not in column_map]
    if missing:
        raise ValueError(
            f"[{label}] column_map is missing required keys {missing}. "
            f"Required: {required}; optional: ['beta', 'rsid']."
        )

    usecols = list(dict.fromkeys(column_map.values()))
    header = pd.read_csv(gwas_path, sep=sep, compression=compression, nrows=0)
    absent = [c for c in usecols if c not in header.columns]
    if absent:
        raise ValueError(
            f"[{label}] mapped columns {absent} are not in {gwas_path}. "
            f"The file's columns are: {header.columns.tolist()}. "
            f"Fix the column mapping and re-run -- do not guess."
        )

    df = pd.read_csv(
        gwas_path,
        sep=sep,
        compression=compression,
        usecols=usecols,
        dtype={column_map["chr"]: str},
    )
    df = df.rename(columns={v: k for k, v in column_map.items()})

    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df = df[df["pval"] < pval_threshold].copy()

    if "beta" not in df.columns:
        df["beta"] = 0.0
        print(
            f"[{label}] no 'beta' column mapped -- filling 0.0 (allele reorientation "
            f"will not flip the effect sign)"
        )

    df["chr"] = df["chr"].astype(str).str.replace(r"^chr", "", regex=True)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    if uppercase_alleles:
        df["ref"] = df["ref"].astype(str).str.upper()
        df["alt"] = df["alt"].astype(str).str.upper()

    df = df.dropna(subset=["chr", "pos", "ref", "alt", "pval"])
    df["pos"] = df["pos"].astype(int)

    if "rsid" not in df.columns:
        df["rsid"] = df["chr"] + ":" + df["pos"].astype(str)

    df["variant"] = (
        df["chr"]
        + ":"
        + df["pos"].astype(str)
        + ":"
        + df["ref"].astype(str)
        + ":"
        + df["alt"].astype(str)
    )
    df = df[CANONICAL]
    print(f"[{label}] {len(df):,} variants below pval<{pval_threshold}")
    return df
