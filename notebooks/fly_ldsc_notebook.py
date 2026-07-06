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


@app.cell
def __(mo):
    mo.md("## 3. Generate cell-type annotations (BED → .annot.gz)")
    return


@app.cell
def __(subprocess, os, all_cell_types, cell_type_beds, python27_path, ldsc27_path, FLY_CHROMS, DGRP_PREFIX):
    os.makedirs("data/annotations", exist_ok=True)
    _env = os.environ.copy()
    _env["PATH"] = f"{ldsc27_path}/bin:" + _env.get("PATH", "")

    for _ct in all_cell_types:
        _all_exist = all(
            os.path.exists(f"data/annotations/{_ct}.{_ch}.annot.gz")
            for _ch in FLY_CHROMS
        )
        if _all_exist:
            print(f"  {_ct}: all annotations exist, skipping")
            continue

        print(f"\nProcessing {_ct}...")
        _bed = cell_type_beds[_ct]
        for _ch in FLY_CHROMS:
            _out = f"data/annotations/{_ct}.{_ch}.annot.gz"
            if os.path.exists(_out):
                continue
            _r = subprocess.run(
                [
                    python27_path, "tools/ldsc/make_annot.py",
                    "--bed-file",   _bed,
                    "--bimfile",    f"{DGRP_PREFIX}.{_ch}.bim",
                    "--annot-file", _out,
                ],
                capture_output=True, text=True, env=_env,
            )
            if _r.returncode != 0:
                print(f"  ERROR chr{_ch}: {_r.stderr[:200]}")
            else:
                print(f"  chr{_ch}", end=" ", flush=True)
        print(f"\n  {_ct} done")

    print("\nAll annotations generated")
    return


@app.cell
def __(mo):
    mo.md("## 4. Calculate LD scores")
    return


@app.cell
def __(subprocess, os, all_cell_types, python27_path, concurrent, multiprocessing, FLY_CHROMS, DGRP_PREFIX):
    os.makedirs("data/ldscores", exist_ok=True)

    def _calc_ld(args):
        ct, ch = args
        _dir = f"data/ldscores/{ct}"
        os.makedirs(_dir, exist_ok=True)
        _out = f"{_dir}/{ct}.{ch}.l2.ldscore.gz"
        if os.path.exists(_out):
            return f"[{ct}] chr{ch} exists"
        try:
            subprocess.run(
                [
                    python27_path, "tools/ldsc/ldsc.py",
                    "--l2",
                    "--bfile",       f"{DGRP_PREFIX}.{ch}",
                    "--ld-wind-kb",  "1000",
                    "--annot",       f"data/annotations/{ct}.{ch}.annot.gz",
                    "--thin-annot",
                    "--out",         f"{_dir}/{ct}.{ch}",
                ],
                check=True, capture_output=True,
            )
            return f"[{ct}] chr{ch} done"
        except subprocess.CalledProcessError as e:
            return f"ERROR [{ct}] chr{ch}: {e.stderr[:200] if e.stderr else ''}"

    _tasks = [(ct, ch) for ct in all_cell_types for ch in FLY_CHROMS]
    _max_workers = min(multiprocessing.cpu_count() - 1, 8)
    print(f"Running LD score calculation ({_max_workers} workers, {len(_tasks)} tasks)...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=_max_workers) as _ex:
        for _res in _ex.map(_calc_ld, _tasks):
            print(_res)

    print("\nAll LD scores calculated")
    return


@app.cell
def __(mo):
    mo.md("""
    ## 5. Phenotype Preparation

    Download longevity phenotype data from DGRPool Study 1 (Arya et al. 2010)
    and average male and female measurements per DGRP line.

    **Source**: https://dgrpool.epfl.ch/studies/1/
    **Phenotype**: Mean longevity (days) per DGRP line, both sexes averaged
    """)
    return


@app.cell
def __(BASE_DIR, pd):
    import urllib.request, gzip as _gzip

    _SUMMARY_URL = "https://dgrpool.epfl.ch/studies/1/get_file?name=summary.tsv"
    _RAW_FILE    = BASE_DIR / "data" / "gwas" / "study1_longevity_raw.tsv"
    _PHENO_FILE  = BASE_DIR / "data" / "gwas" / "lifespan_both_sex.pheno"
    _PHENO_NAME  = "longevity_both_sex"

    (BASE_DIR / "data" / "gwas").mkdir(parents=True, exist_ok=True)

    if not _RAW_FILE.exists():
        print("Downloading Study 1 longevity summary from DGRPool...")
        urllib.request.urlretrieve(_SUMMARY_URL, _RAW_FILE)
    else:
        print(f"Already downloaded: {_RAW_FILE}")   

    with _gzip.open(_RAW_FILE, "rt") as _f:
        _df = pd.read_csv(_f, sep="\t")

    _avg = _df.groupby("DGRP")["mn_Longevity"].mean().reset_index()
    _avg.columns = ["DGRP", _PHENO_NAME]
    _avg["FID"] = "line"
    _avg["IID"] = _avg["DGRP"].str.replace("DGRP_", "", regex=False).str.lstrip("0")
    _avg = _avg[["FID", "IID", _PHENO_NAME]].dropna()

    _fam = pd.read_csv(
        BASE_DIR / "data" / "reference" / "DGRP.2L.fam", sep=" ", header=None,
        names=["FID", "IID", "f", "m", "s", "p"]
    )
    _overlap = set(_avg["IID"].astype(str)) & set(_fam["IID"].astype(str))
    print(f"Lines averaged: {len(_avg)}  |  DGRP2 overlap: {len(_overlap)}")

    _avg.to_csv(_PHENO_FILE, sep="\t", index=False)
    print(f"Saved: {_PHENO_FILE}")
    _avg


@app.cell
def __(mo):
    mo.md("""
    ## 6. Genotype QC

    Filter DGRP2 genotypes per chromosome arm before running GWAS.

    - `--maf 0.01` — remove SNPs with minor allele frequency < 1% (too rare to test reliably with n=132)
    - `--geno 0.05` — remove SNPs missing in more than 5% of lines (low-quality genotyping)

    Cleaned files saved to `data/gwas/tmp/qc/<chrom>.bed/bim/fam`
    """)
    return


@app.cell
def __(BASE_DIR, FLY_CHROMS, subprocess):
    _QC_DIR  = BASE_DIR / "data" / "gwas" / "tmp" / "qc"
    _REF_DIR = BASE_DIR / "data" / "reference"
    _QC_DIR.mkdir(parents=True, exist_ok=True)

    print("Running per-chromosome QC (MAF >= 0.01, SNP missingness <= 0.05)...")
    for _chrom in FLY_CHROMS:
        _out = _QC_DIR / _chrom
        if _out.with_suffix(".bed").exists():
            print(f"  {_chrom}: already done, skipping")
            continue
        _result = subprocess.run([
            "plink2",
            "--bfile",         str(_REF_DIR / f"DGRP.{_chrom}"),
            "--maf",           "0.01",
            "--geno",          "0.05",
            "--allow-extra-chr",
            "--make-bed",
            "--out",           str(_out),
        ], capture_output=True, text=True)
        if _result.returncode != 0:
            print(f"  ERROR on {_chrom}:\n{_result.stderr[-300:]}")
        else:
            _n = sum(1 for _ in open(f"{_out}.bim"))
            print(f"  {_chrom}: {_n:,} SNPs after QC")


@app.cell
def __(mo):
    mo.md("## 7. Process GWAS summary statistics")
    return


@app.cell
def __(GWAS_FILE, SUMSTATS_FILE, os, pd, python27_path):
    if GWAS_FILE is None or not os.path.isfile(GWAS_FILE):
        print("No GWAS file configured — skipping sumstats processing")
        print("When fly GWAS data is available, set GWAS_INPUT_FILE in cell 2.")
    elif SUMSTATS_FILE and os.path.exists(SUMSTATS_FILE):
        print(f"Sumstats already exists: {SUMSTATS_FILE}")
    else:
        import gzip as _gzip
        os.makedirs("data/gwas", exist_ok=True)

        _lower = GWAS_FILE.lower()
        _is_gz = _lower.endswith(".gz") or _lower.endswith(".bgz")
        _open_fn = _gzip.open if _is_gz else open

        with _open_fn(GWAS_FILE, "rt") as _fh:
            for _line in _fh:
                if not _line.startswith("#"):
                    _peek = _line
                    break
        _sep = "\t" if "\t" in _peek else " "

        _df = pd.read_csv(GWAS_FILE, sep=_sep, compression="gzip" if _is_gz else None, low_memory=False)
        _df.columns = [c.lstrip("#").strip() for c in _df.columns]
        print(f"Columns: {list(_df.columns)}")

        _cl = {c.lower(): c for c in _df.columns}
        _snp_col  = _cl.get("snp") or _cl.get("rsid") or _cl.get("variant_id") or _cl.get("id")
        _a1_col   = _cl.get("a1") or _cl.get("effect_allele") or _cl.get("alt")
        _a2_col   = _cl.get("a2") or _cl.get("other_allele") or _cl.get("ref")
        _beta_col = _cl.get("beta") or _cl.get("logor") or _cl.get("b")
        _se_col   = _cl.get("se") or _cl.get("stderr") or _cl.get("standard_error")
        _p_col    = _cl.get("p") or _cl.get("p_value") or _cl.get("pvalue") or _cl.get("p-value")
        _n_col    = _cl.get("n") or _cl.get("n_samples") or _cl.get("neff")
        _z_col    = _cl.get("z") or _cl.get("zscore")

        _rename = {}
        if _snp_col:  _rename[_snp_col]  = "SNP"
        if _a1_col:   _rename[_a1_col]   = "A1"
        if _a2_col:   _rename[_a2_col]   = "A2"
        if _beta_col: _rename[_beta_col] = "BETA"
        if _se_col:   _rename[_se_col]   = "SE"
        if _p_col:    _rename[_p_col]    = "P"
        if _n_col:    _rename[_n_col]    = "N"
        if _z_col:    _rename[_z_col]    = "Z"
        _df = _df.rename(columns=_rename)

        if "Z" not in _df.columns and "BETA" in _df.columns and "SE" in _df.columns:
            _df["Z"] = pd.to_numeric(_df["BETA"], errors="coerce") / pd.to_numeric(_df["SE"], errors="coerce")

        if "N" not in _df.columns:
            print("WARNING: No sample size column found — set N manually below")
            _df["N"] = 1000

        _keep = [c for c in ["SNP", "A1", "A2", "Z", "N"] if c in _df.columns]
        _out  = _df[_keep].dropna(subset=["SNP", "Z"])
        _out.to_csv(SUMSTATS_FILE, sep="\t", index=False, compression="gzip")
        print(f"Written {len(_out):,} variants to {SUMSTATS_FILE}")
    return


@app.cell
def __(mo):
    mo.md("## 6. Build baseline LD scores (fly genome-wide)")
    return


@app.cell
def __(mo):
    mo.md("""
    **Note on fly baseline model:**

    The human pipeline uses a pre-built baseline LD score model (1000G hg38 baseline v1.2).
    For *Drosophila*, no equivalent pre-built baseline exists. You have two options:

    **Option A — No baseline (CTS only)**
    Run `--h2-cts` without `--ref-ld-chr` baseline. This tests enrichment relative to the
    genome-wide LD score, which is simpler but less well-powered.

    **Option B — Compute a genome-wide DGRP baseline**
    Run `ldsc.py --l2` on all DGRP SNPs (no annotation / `--thin-annot`) to get a
    baseline LD score file. Use that as `--ref-ld-chr` in the CTS step.

    Cell 6b below computes Option B. Skip it if you want Option A.
    """)
    return


@app.cell
def __(subprocess, os, python27_path, FLY_CHROMS, DGRP_PREFIX):
    os.makedirs("data/ldscores/baseline", exist_ok=True)

    for _ch in FLY_CHROMS:
        _out = f"data/ldscores/baseline/baseline.{_ch}.l2.ldscore.gz"
        if os.path.exists(_out):
            print(f"  chr{_ch} baseline exists, skipping")
            continue
        print(f"  Computing baseline chr{_ch}...", end=" ", flush=True)
        _r = subprocess.run(
            [
                python27_path, "tools/ldsc/ldsc.py",
                "--l2",
                "--bfile",      f"{DGRP_PREFIX}.{_ch}",
                "--ld-wind-kb", "1000",
                "--out",        f"data/ldscores/baseline/baseline.{_ch}",
            ],
            capture_output=True,
        )
        if _r.returncode != 0:
            print(f"ERROR: {_r.stderr[:200] if _r.stderr else ''}")
        else:
            print("done")

    print("Baseline LD scores ready")
    return


@app.cell
def __(mo):
    mo.md("## 7. Create CTS reference file")
    return


@app.cell
def __(os, all_cell_types, CTS_FILE, FLY_CHROMS):
    os.makedirs("results", exist_ok=True)

    COMPLETED_CELL_TYPES = []
    for _ct in all_cell_types:
        _complete = all(
            os.path.exists(f"data/ldscores/{_ct}/{_ct}.{_ch}.l2.ldscore.gz")
            for _ch in FLY_CHROMS
        )
        if _complete:
            COMPLETED_CELL_TYPES.append(_ct)
        else:
            _missing = [c for c in FLY_CHROMS
                        if not os.path.exists(f"data/ldscores/{_ct}/{_ct}.{c}.l2.ldscore.gz")]
            print(f"  {_ct}: INCOMPLETE — missing chr {_missing}")

    with open(CTS_FILE, "w") as _f:
        for _ct in COMPLETED_CELL_TYPES:
            _f.write(f"{_ct}\tdata/ldscores/{_ct}/{_ct}.\n")

    print(f"CTS file written: {CTS_FILE}  ({len(COMPLETED_CELL_TYPES)} cell types)")
    return (COMPLETED_CELL_TYPES,)


@app.cell
def __(mo):
    mo.md("## 8. Run LDSC cell-type-specific heritability analysis")
    return


@app.cell
def __(CTS_FILE, SUMSTATS_FILE, RESULTS_PREFIX, os, subprocess, python27_path):
    _baseline_exists = os.path.exists("data/ldscores/baseline/baseline.2L.l2.ldscore.gz")

    if SUMSTATS_FILE is None or not os.path.exists(SUMSTATS_FILE):
        print(f"Skipping — sumstats not found: {SUMSTATS_FILE}")
    elif not os.path.exists(CTS_FILE):
        print(f"Skipping — CTS file not found: {CTS_FILE}")
    else:
        os.makedirs("results", exist_ok=True)
        print("Running LDSC CTS analysis...")

        _cmd = [
            python27_path, "tools/ldsc/ldsc.py",
            "--h2-cts",         SUMSTATS_FILE,
            "--ref-ld-chr-cts", CTS_FILE,
            "--out",            RESULTS_PREFIX,
        ]

        if _baseline_exists:
            _cmd += ["--ref-ld-chr", "data/ldscores/baseline/baseline."]
            print("  Using DGRP baseline LD scores")
        else:
            print("  No baseline LD scores — running without baseline (Option A)")

        _r = subprocess.run(_cmd, capture_output=True, text=True)
        print(_r.stdout[-2000:] if _r.stdout else "")
        if _r.returncode != 0:
            print(f"STDERR: {_r.stderr[-1000:]}")
        else:
            print("LDSC CTS analysis complete")
    return


if __name__ == "__main__":
    app.run()
