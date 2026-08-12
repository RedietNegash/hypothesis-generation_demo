#!/usr/bin/env python3
"""
run_hm_through_pipeline.py
==========================
Runs the already-downloaded harmonised file (.h.tsv.gz) through harmonizer.sh,
then LDSC h2-cts, then compares with the existing "skip harmoniser" results.

Advisor's question:
  Path A: .h.tsv.gz → harmonizer.sh → LDSC   (use harmonised file IN pipeline)
  Path B: .h.tsv.gz → LDSC directly           (skip harmoniser)  ← already done
"""
import os, json, glob, subprocess
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

STUDY_ID   = "GCST90387778"
HM_FILE    = f"/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation/{STUDY_ID}/{STUDY_ID}.h.tsv.gz"
RES_SKIP   = f"/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation/{STUDY_ID}/{STUDY_ID}_ldsc_ot.cell_type_results.txt"
STUDY_DIR  = f"/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation/{STUDY_ID}"
N_SAMPLES  = 40282

LDSC_WD    = "/mnt/hdd_1/rediet/hypothesis-generation-demo/notebooks"
W_HM3      = os.path.join(LDSC_WD, "data/OSF/w_hm3.snplist")
CTS_FILE   = os.path.join(LDSC_WD, "data/PASS_AtrialFibrillation_Nielsen2018_sumstats_cell_types.cts")
BASELINE   = "data/OSF/baseline_v1.2/baseline."
LDSCORES   = "data/ldscores/all_merged_cCREs/all_merged_cCREs."
WEIGHTS    = "data/OSF/weights/weights.hm3_noMHC."

HARMONIZER_SCRIPT  = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/harmonizer.sh"
HARMONIZER_REF_DIR = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/data/gwas_harm_ref"
HARMONIZER_REPO    = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser"
NEXTFLOW_ENV = {
    **os.environ,
    "PATH": "/mnt/hdd_1/rediet/hypothesis-generation-demo:/mnt/hdd_1/rediet/jdk-17/bin:" + os.environ.get("PATH", ""),
    "JAVA_HOME": "/mnt/hdd_1/rediet/jdk-17",
}
AMBIGUOUS = [{"A", "T"}, {"C", "G"}]

print("=" * 60)
print(f"Study: {STUDY_ID}  N={N_SAMPLES:,}")
print("Path A: .h.tsv.gz → harmonizer.sh → LDSC")
print("Path B: .h.tsv.gz → LDSC directly  (already done, reusing)")
print("=" * 60)

def load_hm3():
    return pd.read_csv(W_HM3, sep="\t")[["SNP", "A1", "A2"]]


def to_sumstats(df, n):
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"]   = pd.to_numeric(df["SE"],   errors="coerce")
    df["A1"]   = df["A1"].str.upper()
    df["A2"]   = df["A2"].str.upper()
    df["Z"]    = df["BETA"] / df["SE"]
    df["N"]    = n
    df = df.dropna(subset=["SNP", "Z"])
    hm3 = load_hm3()
    df  = df.merge(hm3, on="SNP", suffixes=("", "_hm3"))
    df  = df[(df["A1"] == df["A1_hm3"]) | (df["A1"] == df["A2_hm3"])]
    df  = df.drop(columns=["A1_hm3", "A2_hm3"])
    df  = df[~df.apply(lambda r: {r["A1"], r["A2"]} in AMBIGUOUS, axis=1)]
    return df[["SNP", "A1", "A2", "Z", "N"]].dropna()


def reformat_for_harmonizer(hm_path, study_dir):
    """Reformat the .h.tsv.gz so harmonizer.sh can ingest it."""
    import gzip as _gzip
    out_dir  = os.path.join(study_dir, "gwas_reformatted_hm")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{STUDY_ID}.h.tsv.gz")

    if os.path.exists(out_path):
        print(f"  Already reformatted: {os.path.basename(out_path)}")
        return out_path

    df = pd.read_csv(hm_path, sep="\t", compression="gzip", low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    print(f"  Columns: {list(df.columns)}")

    cl = {c.lower(): c for c in df.columns}
    snp_col  = cl.get("rsid") or cl.get("rs_id") or cl.get("hm_rsid") or cl.get("snp") or cl.get("variant_id")
    chr_col  = cl.get("chromosome") or cl.get("chr") or cl.get("hm_chrom")
    pos_col  = cl.get("base_pair_location") or cl.get("pos") or cl.get("bp")
    a1_col   = cl.get("effect_allele") or cl.get("hm_effect_allele") or cl.get("a1")
    a2_col   = cl.get("other_allele")  or cl.get("hm_other_allele")  or cl.get("a2")
    beta_col = cl.get("beta") or cl.get("hm_beta") or cl.get("effect")
    se_col   = cl.get("standard_error") or cl.get("se") or cl.get("stderr")
    p_col    = cl.get("p_value") or cl.get("p-value") or cl.get("p")

    rename = {c: s for c, s in [(snp_col,"snp"),(chr_col,"chr"),(pos_col,"pos"),
                                  (a1_col,"a1"),(a2_col,"a2"),(beta_col,"beta"),
                                  (se_col,"se"),(p_col,"p")] if c}
    df = df.rename(columns=rename)
    if "chr" in df.columns:
        df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)

    keep      = [c for c in ["chr","pos","snp","a1","a2","beta","se","p"] if c in df.columns]
    drop_cols = [c for c in ["a1","a2","beta","p"] if c in df.columns]
    df[keep].dropna(subset=drop_cols).to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(df):,} variants → {os.path.basename(out_path)}")
    return out_path


def run_harmonizer(reformatted_path, study_dir):
    ref_dir       = os.path.dirname(reformatted_path)
    expected_out  = os.path.join(ref_dir, STUDY_ID)
    harm_work_dir = os.path.join(study_dir, "harmonized_hm")
    os.makedirs(harm_work_dir, exist_ok=True)

    if os.path.exists(expected_out):
        print(f"  Harmonizer output already exists: {expected_out}")
        return expected_out

    log_file  = os.path.join(harm_work_dir, "harmonizer.log")
    input_abs = os.path.abspath(reformatted_path)
    orig_dir  = os.getcwd()
    os.chdir(harm_work_dir)
    try:
        print(f"  Running harmonizer.sh on {os.path.basename(reformatted_path)} ...")
        with open(log_file, "w") as lf:
            result = subprocess.run([
                "bash", HARMONIZER_SCRIPT,
                "--input",     input_abs,
                "--build",     "GRCh38",
                "--ref",       HARMONIZER_REF_DIR,
                "--code-repo", HARMONIZER_REPO,
                "--threshold", "0.99",
            ], check=False, text=True, env=NEXTFLOW_ENV, stdout=lf, stderr=lf)
        print(f"  harmonizer.sh exit code: {result.returncode}")
        return expected_out if os.path.exists(expected_out) else ref_dir
    except Exception as e:
        print(f"  ERROR: {e}")
        return None
    finally:
        os.chdir(orig_dir)


def harmonized_to_sumstats(harm_out_dir, reformatted_path, n):
    final_dir = os.path.join(harm_out_dir, "final")
    files = []
    if os.path.exists(final_dir):
        files = [f for f in os.listdir(final_dir) if f.endswith(".tsv.gz")]
    if not files:
        # Fix: avoid double extension bug for files already named *.tsv.gz
        # Try standard SSF path first, then fall back to the reformatted file itself
        ssf_path = reformatted_path.replace(".tsv.gz", ".tsv.gz")  # no-op, use as-is
        if reformatted_path.endswith(".tsv.gz"):
            ssf_path = reformatted_path  # already correct format, use directly
        else:
            ssf_path = reformatted_path.replace(".gz", ".tsv.gz")
        print(f"  No final/ dir — SSF fallback: {os.path.basename(ssf_path)}")
        if os.path.exists(ssf_path):
            files     = [os.path.basename(ssf_path)]
            final_dir = os.path.dirname(ssf_path)
        else:
            print(f"  WARNING: SSF file not found: {ssf_path}")
            return None
    harmonized_file = os.path.join(final_dir, files[0])
    print(f"  Reading: {os.path.basename(harmonized_file)}")
    df = pd.read_csv(harmonized_file, sep="\t", compression="gzip")
    df = df.rename(columns={"rsid": "SNP", "effect_allele": "A1", "other_allele": "A2",
                             "beta": "BETA", "standard_error": "SE", "p_value": "P"})
    if "N" not in df.columns:
        df["N"] = n
    return to_sumstats(df, n)


def get_python27():
    j    = subprocess.run(["conda", "env", "list", "--json"], capture_output=True, text=True, check=True)
    envs = json.loads(j.stdout)["envs"]
    return os.path.join([e for e in envs if "ldsc27" in e][0], "bin", "python")


def run_ldsc(sumstats_path, out_prefix, python27):
    results_file = f"{out_prefix}.cell_type_results.txt"
    if os.path.exists(results_file):
        print(f"  Results already exist: {os.path.basename(results_file)}")
        return results_file
    cmd = [
        python27, os.path.join(LDSC_WD, "tools/ldsc/ldsc.py"),
        "--h2-cts",         sumstats_path,
        "--ref-ld-chr",     f"{BASELINE},{LDSCORES}",
        "--ref-ld-chr-cts", CTS_FILE,
        "--w-ld-chr",       os.path.join(LDSC_WD, WEIGHTS),
        "--out",            out_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=LDSC_WD)
    if result.returncode != 0:
        print(f"  LDSC ERROR:\n{result.stderr[-500:]}")
        return None
    return results_file



python27 = get_python27()
print(f"Python 2.7: {python27}\n")


print("── Path A: .h.tsv.gz → harmonizer.sh → LDSC ──")
reformatted = reformat_for_harmonizer(HM_FILE, STUDY_DIR)
harm_out    = run_harmonizer(reformatted, STUDY_DIR)
ss_a        = harmonized_to_sumstats(harm_out, reformatted, N_SAMPLES)

if ss_a is None or len(ss_a) < 1000:
    print(f"Too few SNPs ({len(ss_a) if ss_a else 0}) — aborting.")
    exit(1)

ss_a_path = os.path.join(STUDY_DIR, f"{STUDY_ID}_hm_via_harmonizer.sumstats.gz")
ss_a.to_csv(ss_a_path, sep="\t", index=False, compression="gzip")
print(f"  Written {len(ss_a):,} SNPs → {os.path.basename(ss_a_path)}\n")

print("── Running LDSC h2-cts (Path A) ──")
res_a = run_ldsc(ss_a_path, os.path.join(STUDY_DIR, f"{STUDY_ID}_ldsc_hm_via_harmonizer"), python27)


print(f"\n── Path B: reusing existing results ({os.path.basename(RES_SKIP)}) ──")
res_b = RES_SKIP

if res_a and os.path.exists(res_b):
    r_a = pd.read_csv(res_a, sep="\t").sort_values("Coefficient_P_value")
    r_b = pd.read_csv(res_b, sep="\t").sort_values("Coefficient_P_value")
    _, r_a["FDR"] = fdrcorrection(r_a["Coefficient_P_value"].fillna(1))
    _, r_b["FDR"] = fdrcorrection(r_b["Coefficient_P_value"].fillna(1))

    merged = r_a[["Name", "Coefficient_P_value"]].merge(
        r_b[["Name", "Coefficient_P_value"]], on="Name", suffixes=("_a", "_b")
    )
    rho, p = stats.spearmanr(merged["Coefficient_P_value_a"], merged["Coefficient_P_value_b"])

    print(f"\n{'='*60}")
    print(f"  RESULT — {STUDY_ID}")
    print(f"{'='*60}")
    print(f"  Spearman ρ (222 cell types): {rho:.4f}  (p={p:.2e})")
    print(f"  Sig p<0.05  Path A (via harmoniser): {(r_a['Coefficient_P_value']<0.05).sum():>3}")
    print(f"  Sig p<0.05  Path B (skip harmoniser): {(r_b['Coefficient_P_value']<0.05).sum():>3}")
    print(f"\n  Top 5 — Path A (.h.tsv.gz → harmonizer.sh):")
    print(r_a[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))
    print(f"\n  Top 5 — Path B (.h.tsv.gz → LDSC directly):")
    print(r_b[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))

    if rho > 0.95:
        verdict = "CONSISTENT — can skip harmonizer.sh when OT file is available"
    elif rho > 0.80:
        verdict = " MOSTLY CONSISTENT — minor differences"
    else:
        verdict = "DISCORDANT — harmonizer.sh changes results even on pre-harmonised file"
    print(f"\n  Verdict: {verdict}")
else:
    print("ERROR: missing results file(s)")

print("\nDone.")
