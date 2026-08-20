#!/usr/bin/env python3
import gzip
import re
import shutil
import subprocess
from pathlib import Path

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
FLY_CHROMS = ["2L", "2R", "3L", "3R", "4", "X"]
CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "X", "23": "X"}

MAGMA_BIN = BASE_DIR / "tools" / "magma" / "magma"

GTF_FILE = BASE_DIR / "data" / "genes" / "Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
GTF_SHARED_PATH = Path(
    "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/dmel/gencode"
    "/Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
)
GTF_URL = (
    "https://ftp.ensembl.org/pub/release-112/gtf/drosophila_melanogaster/"
    "Drosophila_melanogaster.BDGP6.46.112.gtf.gz"
)

MAGMA_DIR = BASE_DIR / "data" / "magma"
GENE_LOC_FILE = MAGMA_DIR / "gene_loc.txt"
SNP_LOC_FILE = MAGMA_DIR / "snp_loc.txt"
ANNOT_PREFIX = MAGMA_DIR / "dgrp_lifespan_female"
ANNOT_FILE = MAGMA_DIR / "dgrp_lifespan_female.genes.annot"

MERGED_QC_SOURCE = BASE_DIR / "data" / "gwas" / "tmp" / "merged_qc"
MAGMA_BFILE = MAGMA_DIR / "merged_qc"
PHENO_FILE = BASE_DIR / "data" / "gwas" / "lifespan_female.pheno"
PHENO_NAME = "S18_1537_F"
GENE_OUT_PREFIX = MAGMA_DIR / "dgrp_lifespan_female_gene"
GENE_OUT_FILE = MAGMA_DIR / "dgrp_lifespan_female_gene.genes.out"
GENE_RAW_FILE = MAGMA_DIR / "dgrp_lifespan_female_gene.genes.raw"

PEAKS_DIR = BASE_DIR / "data" / "peaks"
GENE_SETS_FILE = MAGMA_DIR / "gene_sets.txt"
GENESET_OUT_PREFIX = MAGMA_DIR / "dgrp_lifespan_female_geneset"
GENESET_OUT_FILE = MAGMA_DIR / "dgrp_lifespan_female_geneset.gsa.out"

WINDOW_KB = "5"
ANNOT_PREFIX_W5 = MAGMA_DIR / "dgrp_lifespan_female_w5"
ANNOT_FILE_W5 = MAGMA_DIR / "dgrp_lifespan_female_w5.genes.annot"
GENE_OUT_PREFIX_W5 = MAGMA_DIR / "dgrp_lifespan_female_gene_w5"
GENE_OUT_FILE_W5 = MAGMA_DIR / "dgrp_lifespan_female_gene_w5.genes.out"
GENE_RAW_FILE_W5 = MAGMA_DIR / "dgrp_lifespan_female_gene_w5.genes.raw"
GENESET_OUT_PREFIX_W5 = MAGMA_DIR / "dgrp_lifespan_female_geneset_w5"
GENESET_OUT_FILE_W5 = MAGMA_DIR / "dgrp_lifespan_female_geneset_w5.gsa.out"


def ensure_gtf():
    if GTF_FILE.exists():
        return
    GTF_FILE.parent.mkdir(parents=True, exist_ok=True)
    if GTF_SHARED_PATH.exists():
        shutil.copy(GTF_SHARED_PATH, GTF_FILE)
    else:
        subprocess.run(["curl", "-sL", "-o", str(GTF_FILE), GTF_URL], check=True)


def build_gene_location_file():
    if GENE_LOC_FILE.exists():
        print(f"Gene location file already exists: {GENE_LOC_FILE}")
        return

    ensure_gtf()
    print("Parsing genome-wide gene coordinates from GTF ...", flush=True)
    MAGMA_DIR.mkdir(parents=True, exist_ok=True)
    keep = set(FLY_CHROMS)
    name_re = re.compile(r'gene_name "([^"]+)"')
    id_re = re.compile(r'gene_id "([^"]+)"')

    n = 0
    with gzip.open(GTF_FILE, "rt") as fh, open(GENE_LOC_FILE, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene" or f[0] not in keep:
                continue
            m = name_re.search(f[8]) or id_re.search(f[8])
            if not m:
                continue
            gene, chrom, start, stop, strand = m.group(1), f[0], f[3], f[4], f[6]
            out.write(f"{gene}\t{CHROM_MAP[chrom]}\t{start}\t{stop}\t{strand}\n")
            n += 1

    print(f"  {n:,} genes written to {GENE_LOC_FILE}")


def build_snp_location_file():
    if SNP_LOC_FILE.exists():
        print(f"SNP location file already exists: {SNP_LOC_FILE}")
        return

    MAGMA_DIR.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(SNP_LOC_FILE, "w") as out:
        for ch in FLY_CHROMS:
            bim_file = BASE_DIR / "data" / "reference" / f"DGRP.{ch}.bim"
            with open(bim_file) as fh:
                for line in fh:
                    chrom, snp, _cm, pos, _a1, _a2 = line.rstrip("\n").split("\t")
                    out.write(f"{snp}\t{CHROM_MAP[chrom]}\t{pos}\n")
                    n += 1

    print(f"  {n:,} SNPs written to {SNP_LOC_FILE}")


def run_annotate():
    if ANNOT_FILE.exists():
        print(f"Annotation file already exists: {ANNOT_FILE}")
        return

    subprocess.run([
        str(MAGMA_BIN),
        "--annotate", "nonhuman",
        "--snp-loc", str(SNP_LOC_FILE),
        "--gene-loc", str(GENE_LOC_FILE),
        "--out", str(ANNOT_PREFIX),
    ], check=True)


def prepare_magma_bfile():
    bim_out = MAGMA_BFILE.with_suffix(".bim")
    if bim_out.exists():
        print(f"MAGMA bfile already prepared: {MAGMA_BFILE}")
        return

    MAGMA_DIR.mkdir(parents=True, exist_ok=True)
    for ext in [".bed", ".fam"]:
        dst = MAGMA_BFILE.with_suffix(ext)
        if not dst.exists():
            dst.symlink_to(MERGED_QC_SOURCE.with_suffix(ext).resolve())

    n = 0
    with open(MERGED_QC_SOURCE.with_suffix(".bim")) as fh, open(bim_out, "w") as out:
        for line in fh:
            chrom, snp, cm, pos, a1, a2 = line.rstrip("\n").split("\t")
            out.write(f"{CHROM_MAP[chrom]}\t{snp}\t{cm}\t{pos}\t{a1}\t{a2}\n")
            n += 1

    print(f"  {n:,} SNPs remapped into {bim_out}")


def run_gene_analysis():
    if GENE_OUT_FILE.exists():
        print(f"Gene analysis results already exist: {GENE_OUT_FILE}")
        return

    prepare_magma_bfile()
    subprocess.run([
        str(MAGMA_BIN),
        "--bfile", str(MAGMA_BFILE),
        "--pheno", f"file={PHENO_FILE}", f"use={PHENO_NAME}",
        "--covar", "chrX-use-sex=0",
        "--gene-annot", str(ANNOT_FILE),
        "--out", str(GENE_OUT_PREFIX),
    ], check=True)


def build_gene_sets_file():
    if GENE_SETS_FILE.exists():
        print(f"Gene sets file already exists: {GENE_SETS_FILE}")
        return

    MAGMA_DIR.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(GENE_SETS_FILE, "w") as out:
        for bed in sorted(PEAKS_DIR.glob("*.bed")):
            cell_type = bed.stem
            genes = []
            seen = set()
            with open(bed) as fh:
                for line in fh:
                    gene = line.rstrip("\n").split("\t")[3]
                    if gene not in seen:
                        seen.add(gene)
                        genes.append(gene)
            out.write(cell_type + "\t" + "\t".join(genes) + "\n")
            n += 1

    print(f"  {n:,} cell-type gene sets written to {GENE_SETS_FILE}")


def run_geneset_analysis():
    if GENESET_OUT_FILE.exists():
        print(f"Gene-set results already exist: {GENESET_OUT_FILE}")
        return

    build_gene_sets_file()
    subprocess.run([
        str(MAGMA_BIN),
        "--gene-results", str(GENE_RAW_FILE),
        "--set-annot", str(GENE_SETS_FILE),
        "--out", str(GENESET_OUT_PREFIX),
    ], check=True)


def run_annotate_window():
    if ANNOT_FILE_W5.exists():
        print(f"Windowed annotation file already exists: {ANNOT_FILE_W5}")
        return

    subprocess.run([
        str(MAGMA_BIN),
        "--annotate", f"window={WINDOW_KB}", "nonhuman",
        "--snp-loc", str(SNP_LOC_FILE),
        "--gene-loc", str(GENE_LOC_FILE),
        "--out", str(ANNOT_PREFIX_W5),
    ], check=True)


def run_gene_analysis_window():
    if GENE_OUT_FILE_W5.exists():
        print(f"Windowed gene analysis results already exist: {GENE_OUT_FILE_W5}")
        return

    prepare_magma_bfile()
    subprocess.run([
        str(MAGMA_BIN),
        "--bfile", str(MAGMA_BFILE),
        "--pheno", f"file={PHENO_FILE}", f"use={PHENO_NAME}",
        "--covar", "chrX-use-sex=0",
        "--gene-annot", str(ANNOT_FILE_W5),
        "--out", str(GENE_OUT_PREFIX_W5),
    ], check=True)


def run_geneset_analysis_window():
    if GENESET_OUT_FILE_W5.exists():
        print(f"Windowed gene-set results already exist: {GENESET_OUT_FILE_W5}")
        return

    build_gene_sets_file()
    subprocess.run([
        str(MAGMA_BIN),
        "--gene-results", str(GENE_RAW_FILE_W5),
        "--set-annot", str(GENE_SETS_FILE),
        "--out", str(GENESET_OUT_PREFIX_W5),
    ], check=True)


if __name__ == "__main__":
    build_gene_location_file()
    build_snp_location_file()
    run_annotate()
    run_gene_analysis()
    run_geneset_analysis()
    run_annotate_window()
    run_gene_analysis_window()
    run_geneset_analysis_window()
