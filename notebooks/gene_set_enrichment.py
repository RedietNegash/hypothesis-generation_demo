#!/usr/bin/env python3
import gzip
import re
import shutil
import subprocess
from pathlib import Path

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
FLY_CHROMS = ["2L", "2R", "3L", "3R", "4", "X"]

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
            out.write(f"{gene}\t{chrom}\t{start}\t{stop}\t{strand}\n")
            n += 1

    print(f"  {n:,} genes written to {GENE_LOC_FILE}")


def build_snp_location_file():
    raise NotImplementedError


def run_annotate():
    raise NotImplementedError


def run_gene_analysis():
    raise NotImplementedError


def run_geneset_analysis():
    raise NotImplementedError


if __name__ == "__main__":
    build_gene_location_file()
