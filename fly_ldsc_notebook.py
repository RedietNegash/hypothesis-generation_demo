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


if __name__ == "__main__":
    app.run()
