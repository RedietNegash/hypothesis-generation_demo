# -*- coding: utf-8 -*-
import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import os
    import re
    import subprocess
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import glob
    import concurrent.futures
    import multiprocessing
    return mo, os, re, subprocess, pd, np, Path, glob, concurrent, multiprocessing


@app.cell
def __(mo):
    mo.md("""
    # Fly (Drosophila) LDSC Cell-Type-Specific Heritability Analysis

    This notebook runs LDSC partitioned heritability (cell-type-specific) on
    *Drosophila melanogaster* data using the DGRP (Drosophila Genetic Reference Panel)
    as the LD reference panel.

    **Genome assembly**: dm6

    **Chromosomes used**: 2L, 2R, 3L, 3R, 4, X

    **Reference panel**: DGRP2 plink files (one per chromosome) in `data/reference/`

    **Cell-type peaks**: one BED file per cell type in `data/peaks/`

    Set `GWAS_INPUT_FILE` below when fly GWAS sumstats are available.
    """)
    return


@app.cell
def __(mo, os):
    BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")

    GWAS_INPUT_FILE = mo.ui.text(
        value="data/gwas/",
        label="GWAS input file path (leave empty to skip GWAS steps)",
        full_width=True,
    )

    FLY_CHROMS = ["2L", "2R", "3L", "3R", "4", "X"]

    DGRP_PREFIX = "data/reference/DGRP"

    mo.vstack([
        mo.md("### Configuration"),
        GWAS_INPUT_FILE,
        mo.md(f"Base dir   : `{BASE_DIR}`"),
        mo.md(f"Chromosomes: `{', '.join(FLY_CHROMS)}`"),
        mo.md(f"DGRP prefix: `{DGRP_PREFIX}.<chrom>.bed/bim/fam`"),
    ])
    return BASE_DIR, GWAS_INPUT_FILE, FLY_CHROMS, DGRP_PREFIX


@app.cell
def __(GWAS_INPUT_FILE, os, re):
    _path = GWAS_INPUT_FILE.value.strip()
    if _path and os.path.isfile(_path):
        _basename = os.path.basename(_path)
        _no_ext = _basename
        for _ext in [".tsv.gz", ".txt.gz", ".gz", ".tsv", ".txt", ".csv", ".bgz"]:
            if _no_ext.endswith(_ext):
                _no_ext = _no_ext[: -len(_ext)]
                break
        GWAS_STEM      = re.sub(r"[^A-Za-z0-9_\-]", "_", _no_ext)
        GWAS_FILE      = _path
        SUMSTATS_FILE  = f"data/gwas/{GWAS_STEM}.sumstats.gz"
        CTS_FILE       = f"data/{GWAS_STEM}_cell_types.cts"
        RESULTS_PREFIX = f"results/{GWAS_STEM}_CellTypeSpecific"
        print(f"GWAS stem      : {GWAS_STEM}")
        print(f"GWAS file      : {GWAS_FILE}")
        print(f"Sumstats file  : {SUMSTATS_FILE}")
        print(f"Results prefix : {RESULTS_PREFIX}")
    else:
        GWAS_STEM      = None
        GWAS_FILE      = None
        SUMSTATS_FILE  = None
        CTS_FILE       = "data/cell_types.cts"
        RESULTS_PREFIX = "results/CellTypeSpecific"
        print("No GWAS file specified — GWAS steps will be skipped")
    return GWAS_STEM, GWAS_FILE, SUMSTATS_FILE, CTS_FILE, RESULTS_PREFIX


@app.cell
def __(mo):
    mo.md("## 0. Setup: download and configure LDSC")
    return


@app.cell
def __(Path, subprocess, os):
    import json as _json

    TOOLS_DIR = Path("tools")
    LDSC_DIR  = TOOLS_DIR / "ldsc"
    TOOLS_DIR.mkdir(exist_ok=True)

    _env_check = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if "ldsc27" not in _env_check.stdout:
        subprocess.run(["conda", "create", "-n", "ldsc27", "python=2.7", "-y"], check=True)

    _conda_json = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    _envs = _json.loads(_conda_json.stdout)["envs"]
    ldsc27_path = [e for e in _envs if "ldsc27" in e][0]

    if not os.path.exists(os.path.join(ldsc27_path, "bin", "bedtools")):
        subprocess.run(
            ["conda", "install", "-n", "ldsc27", "-c", "bioconda", "bedtools", "-y"],
            check=True,
        )

    if not LDSC_DIR.exists():
        subprocess.run(
            ["git", "clone", "https://github.com/bulik/ldsc.git", str(LDSC_DIR)],
            check=True,
        )

    _check_np = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "python"), "-c", "import numpy"],
        capture_output=True,
    )
    if _check_np.returncode != 0:
        subprocess.run(
            ["conda", "install", "-n", "ldsc27", "-y", "openssl=1.0.2", "-c", "conda-forge"],
            check=True,
        )
        subprocess.run(
            ["conda", "install", "-n", "ldsc27", "-y",
             "numpy", "scipy", "pandas", "bitarray", "-c", "conda-forge"],
            check=True,
        )
        subprocess.run(
            ["conda", "install", "-n", "ldsc27", "-y",
             "pybedtools", "pysam=0.15.3", "-c", "bioconda", "-c", "conda-forge"],
            check=True,
        )

    subprocess.run(["chmod", "+x", str(LDSC_DIR / "ldsc.py")], check=True)
    subprocess.run(["chmod", "+x", str(LDSC_DIR / "make_annot.py")], check=True)
    python27_path = os.path.join(ldsc27_path, "bin", "python")

    print(f"LDSC environment ready")
    print(f"Python 2.7 : {python27_path}")
    print(f"LDSC dir   : {LDSC_DIR}")
    return TOOLS_DIR, LDSC_DIR, ldsc27_path, python27_path


@app.cell
def __(mo):
    mo.md("## 1. Validate DGRP reference files")
    return


@app.cell
def __(os, FLY_CHROMS, DGRP_PREFIX):
    os.makedirs("data/reference", exist_ok=True)

    _missing = []
    for _ch in FLY_CHROMS:
        for _ext in [".bed", ".bim", ".fam"]:
            _f = f"{DGRP_PREFIX}.{_ch}{_ext}"
            if not os.path.exists(_f):
                _missing.append(_f)

    if _missing:
        print("WARNING: The following DGRP reference files are missing:")
        for _f in _missing:
            print(f"  {_f}")
        print("\nTo generate DGRP plink files:")
        print("  1. Download DGRP2 VCF from http://dgrp2.gnets.ncsu.edu/")
        print("  2. Split by chromosome: bcftools view -r {chrom} dgrp2.vcf.gz | ...")
        print("  3. Convert to plink: plink --vcf dgrp2.{chrom}.vcf.gz --make-bed --out data/reference/DGRP.{chrom}")
        dgrp_ready = False
    else:
        print(f"All DGRP reference files present for chromosomes: {FLY_CHROMS}")
        for _ch in FLY_CHROMS:
            _bim = f"{DGRP_PREFIX}.{_ch}.bim"
            _n = sum(1 for _ in open(_bim))
            print(f"  chr{_ch}: {_n:,} SNPs")
        dgrp_ready = True

    return (dgrp_ready,)


@app.cell
def __(mo):
    mo.md("## 2. Discover cell-type BED files")
    return


@app.cell
def __(os, re):
    os.makedirs("data/peaks", exist_ok=True)

    def _sanitize(name):
        return re.sub(r"[^A-Za-z0-9_\-]", "_", name)

    _bed_files = [
        f for f in os.listdir("data/peaks")
        if f.endswith(".bed") and not f.startswith(".")
    ]

    cell_type_beds = {}
    _seen = {}
    for _f in sorted(_bed_files):
        _raw  = os.path.splitext(_f)[0]
        _safe = _sanitize(_raw)
        if _safe in _seen:
            _safe += "_2"
        _seen[_safe] = _raw
        cell_type_beds[_safe] = os.path.join("data/peaks", _f)

    all_cell_types = sorted(cell_type_beds.keys())

    if not all_cell_types:
        print("WARNING: No BED files found in data/peaks/")
        print("  Place one BED file per cell type: data/peaks/<cell_type>.bed")
        print("  BED format: chrom (e.g. chr2L), start, end (dm6 coordinates)")
    else:
        print(f"{len(all_cell_types)} cell types found:")
        for _ct in all_cell_types:
            _n = sum(1 for _ in open(cell_type_beds[_ct]))
            print(f"  {_ct}: {_n:,} peaks")

    return (all_cell_types, cell_type_beds)


if __name__ == "__main__":
    app.run()
