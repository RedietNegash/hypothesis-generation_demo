#!/usr/bin/env python3
import gzip
import re
from pathlib import Path

import pandas as pd

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")

GTF_FILE = BASE_DIR / "data" / "genes" / "Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
SUSIE_DIR = BASE_DIR / "data" / "finemap" / "female" / "susie_results"
OUT_FILE = BASE_DIR / "data" / "finemap" / "female" / "female_finemap_gene_mapping.tsv"

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
GENE_BIOTYPE_RE = re.compile(r'gene_biotype "([^"]+)"')


def extract_attr(pattern: re.Pattern, attrs: str, default: str = "") -> str:
    m = pattern.search(attrs)
    return m.group(1) if m else default


def load_genes() -> pd.DataFrame:
    rows = []
    with gzip.open(GTF_FILE, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            attrs = fields[8]
            gene_id = extract_attr(GENE_ID_RE, attrs)
            rows.append({
                "CHR": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "gene_id": gene_id,
                "gene_name": extract_attr(GENE_NAME_RE, attrs, default=gene_id),
                "gene_biotype": extract_attr(GENE_BIOTYPE_RE, attrs),
            })
    return pd.DataFrame(rows)


def nearest_genes(genes: pd.DataFrame, chrom: str, pos: int) -> pd.DataFrame:
    on_chrom = genes[genes["CHR"] == chrom]

    overlapping = on_chrom[(on_chrom["start"] <= pos) & (on_chrom["end"] >= pos)]
    if len(overlapping) > 0:
        out = overlapping.copy()
        out["distance_bp"] = 0
        out["relation"] = "overlapping"
        return out

    gap_after_gene = pos - on_chrom["end"]
    gap_before_gene = on_chrom["start"] - pos
    distance = pd.concat([gap_after_gene, gap_before_gene], axis=1).max(axis=1)
    nearest_idx = distance.idxmin()

    out = on_chrom.loc[[nearest_idx]].copy()
    out["distance_bp"] = int(distance.loc[nearest_idx])
    out["relation"] = "nearest"
    return out


def load_finemap_variants() -> pd.DataFrame:
    frames = []
    for f in sorted(SUSIE_DIR.glob("*_susie.tsv")):
        df = pd.read_csv(f, sep="\t")
        df["locus"] = f.stem.replace("_susie", "")
        frames.append(df)
    variants = pd.concat(frames, ignore_index=True)

    in_credible_set = variants["CS"].notna()
    top_per_locus = variants.loc[variants.groupby("locus")["PIP"].idxmax()]
    top_per_locus = top_per_locus[~top_per_locus["SNP"].isin(variants.loc[in_credible_set, "SNP"])]

    selected = pd.concat([variants[in_credible_set], top_per_locus])
    selected["variant_role"] = ["credible_set"] * in_credible_set.sum() + ["top_pip_no_cs"] * len(top_per_locus)
    return selected.sort_values(["locus", "PIP"], ascending=[True, False])


if __name__ == "__main__":
    genes = load_genes()
    variants = load_finemap_variants()

    rows = []
    for _, v in variants.iterrows():
        chrom = v["SNP"].split("_")[0]
        matches = nearest_genes(genes, chrom, int(v["POS"]))
        for _, g in matches.iterrows():
            rows.append({
                "locus": v["locus"],
                "SNP": v["SNP"],
                "CHR": chrom,
                "POS": v["POS"],
                "PIP": v["PIP"],
                "CS": v["CS"],
                "variant_role": v["variant_role"],
                "gene_id": g["gene_id"],
                "gene_name": g["gene_name"],
                "gene_biotype": g["gene_biotype"],
                "relation": g["relation"],
                "distance_bp": g["distance_bp"],
            })

    out = pd.DataFrame(rows)
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_FILE, sep="\t", index=False)

    print(out.to_string(index=False))
    print(f"\nSaved: {OUT_FILE}")
