#!/usr/bin/env python3
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
PEAKS_DIR = BASE_DIR / "data" / "peaks"
REF_DIR = BASE_DIR / "data" / "reference"
ANNOT_DIR = BASE_DIR / "data" / "annotations"
LDSC = BASE_DIR / "tools" / "ldsc" / "make_annot.py"

CHROMS = ["2L", "2R", "3L", "3R", "4", "X"]
PYTHON = "conda run -n ldsc_py3 python"
MAX_WORKERS = 8

ANNOT_DIR.mkdir(parents=True, exist_ok=True)

bed_files = sorted(PEAKS_DIR.glob("*.bed"))
cell_types = [f.stem for f in bed_files]

tasks = []
for bed in bed_files:
    ct = bed.stem
    for chrom in CHROMS:
        out = ANNOT_DIR / f"{ct}.{chrom}.annot.gz"
        if out.exists():
            continue
        tasks.append((ct, chrom, bed, out))

print(f"{len(cell_types)} cell types, {len(tasks)} annotation jobs to run ({MAX_WORKERS} workers)")

def run_job(args):
    ct, chrom, bed, out = args
    bim = REF_DIR / f"DGRP.{chrom}.bim"
    cmd = (
        f"{PYTHON} {LDSC} "
        f"--bed-file {bed} "
        f"--bimfile {bim} "
        f"--annot-file {out}"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return ct, chrom, False, result.stderr[-200:]
    return ct, chrom, True, None

done = 0
failed = []
with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
    futures = {ex.submit(run_job, t): t for t in tasks}
    for fut in as_completed(futures):
        ct, chrom, ok, err = fut.result()
        done += 1
        if ok:
            print(f"[{done}/{len(tasks)}] {ct}.{chrom} done", flush=True)
        else:
            failed.append((ct, chrom, err))
            print(f"[{done}/{len(tasks)}] ERROR {ct}.{chrom}: {err}", flush=True)

print(f"\nFinished. {len(tasks) - len(failed)}/{len(tasks)} succeeded.")
if failed:
    print(f"{len(failed)} failed jobs:")
    for ct, chrom, err in failed:
        print(f"  {ct}.{chrom}: {err}")
