#!/usr/bin/env python3
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
EXPR_FILE = Path(
    "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/afca"
    "/afca_afca_annotation_group_by_mean.tsv.gz"
)
GTF_FILE = Path(
    "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/dmel/gencode"
    "/Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
)
OUT_DIR = BASE_DIR / "data" / "peaks"

TOP_N = 500
KEEP_CHROMS = {"2L", "2R", "3L", "3R", "4", "X"}

print("Parsing GTF ...", flush=True)

gene_name_re = re.compile(r'gene_name "([^"]+)"')

gtf_rows = []
with gzip.open(GTF_FILE, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[2] != "gene":
            continue
        chrom = fields[0]
        if chrom not in KEEP_CHROMS:
            continue
        start = int(fields[3]) - 1
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]
        m = gene_name_re.search(attrs)
        if not m:
            continue
        gene_name = m.group(1)
        gtf_rows.append((gene_name, chrom, start, end, strand))

gtf = (
    pd.DataFrame(gtf_rows, columns=["gene_name", "chrom", "start", "end", "strand"])
    .drop_duplicates("gene_name")
    .set_index("gene_name")
)
print(f"  {len(gtf):,} unique genes parsed from GTF", flush=True)

print("Reading expression matrix ...", flush=True)
expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0, compression="gzip")
print(f"  Shape: {expr.shape} (genes × samples)", flush=True)

expr.columns = pd.Index(expr.columns, name="sample")
cell_type_map = {col: re.sub(r"_\d+$", "", col) for col in expr.columns}
cell_types_order = list(dict.fromkeys(cell_type_map.values()))

print(f"  {len(cell_types_order)} unique cell types", flush=True)

expr_by_ct = expr.T.copy()
expr_by_ct.index = pd.Index([cell_type_map[s] for s in expr_by_ct.index])
expr_mean = expr_by_ct.groupby(level=0).mean().T

OUT_DIR.mkdir(parents=True, exist_ok=True)

missing_genes: set[str] = set()
written = 0

for ct in cell_types_order:
    series = expr_mean[ct].dropna()
    top_genes = series.nlargest(TOP_N).index.tolist()

    rows = []
    for g in top_genes:
        if g not in gtf.index:
            missing_genes.add(g)
            continue
        rec = gtf.loc[g]
        rows.append((rec["chrom"], rec["start"], rec["end"], g, 0, rec["strand"]))

    if not rows:
        print(f"  WARNING: no GTF hits for {ct!r}, skipping", flush=True)
        continue

    bed = pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "score", "strand"])
    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)

    safe_name = re.sub(r"[^A-Za-z0-9_\-]", "_", ct).strip("_")
    out_path = OUT_DIR / f"{safe_name}.bed"
    bed.to_csv(out_path, sep="\t", header=False, index=False)
    written += 1

print(f"\nDone. Wrote {written} BED files to {OUT_DIR}/")
if missing_genes:
    print(
        f"  {len(missing_genes)} gene symbols had no GTF match "
        f"(e.g. {sorted(missing_genes)[:5]})"
    )
