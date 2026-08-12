#!/usr/bin/env python3
"""
run_raw_direct.py
=================
Inverse question — Group B:
  Path A: raw file → LDSC directly   (no harmonisation)   ← THIS SCRIPT
  Path B: raw file → harmonizer.sh → LDSC                 ← already done (_ldsc_harmonized)

If ρ < 0.95 → harmonisation changes results on raw files
             (strengthens case for using OT pre-harmonised files)
"""
import os, json, subprocess
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

STUDY_ID  = "GCST90387778"
N_SAMPLES = 40282
STUDY_DIR = f"/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation/{STUDY_ID}"
RAW_FILE  = os.path.join(STUDY_DIR, f"{STUDY_ID}_raw_original.gz")
RES_HM    = os.path.join(STUDY_DIR, f"{STUDY_ID}_ldsc_harmonized.cell_type_results.txt")

LDSC_WD  = "/mnt/hdd_1/rediet/hypothesis-generation-demo/notebooks"
W_HM3    = os.path.join(LDSC_WD, "data/OSF/w_hm3.snplist")
CTS_FILE = os.path.join(LDSC_WD, "data/PASS_AtrialFibrillation_Nielsen2018_sumstats_cell_types.cts")
BASELINE = "data/OSF/baseline_v1.2/baseline."
LDSCORES = "data/ldscores/all_merged_cCREs/all_merged_cCREs."
WEIGHTS  = "data/OSF/weights/weights.hm3_noMHC."

AMBIGUOUS = [{"A","T"}, {"C","G"}]

print("=" * 60)
print(f"Study: {STUDY_ID}  N={N_SAMPLES:,}")
print("Inverse question (Group B):")
print("  Path A: raw file → LDSC directly   ← running now")
print("  Path B: raw file → harmonizer.sh → LDSC  (reusing _ldsc_harmonized)")
print("=" * 60)


def load_hm3():
    return pd.read_csv(W_HM3, sep="\t")[["SNP","A1","A2"]]


def to_sumstats(df, n):
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"]   = pd.to_numeric(df["SE"],   errors="coerce")
    df["A1"]   = df["A1"].str.upper()
    df["A2"]   = df["A2"].str.upper()
    df["Z"]    = df["BETA"] / df["SE"]
    df["N"]    = n
    df = df.dropna(subset=["SNP","Z"])
    hm3 = load_hm3()
    df  = df.merge(hm3, on="SNP", suffixes=("","_hm3"))
    df  = df[(df["A1"]==df["A1_hm3"]) | (df["A1"]==df["A2_hm3"])]
    df  = df.drop(columns=["A1_hm3","A2_hm3"])
    df  = df[~df.apply(lambda r: {r["A1"],r["A2"]} in AMBIGUOUS, axis=1)]
    return df[["SNP","A1","A2","Z","N"]].dropna()


def get_python27():
    j    = subprocess.run(["conda","env","list","--json"], capture_output=True, text=True, check=True)
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


# ── Path A: raw → LDSC directly ──
print("\n── Path A: raw → LDSC directly ──")
ss_path = os.path.join(STUDY_DIR, f"{STUDY_ID}_raw_direct.sumstats.gz")

if os.path.exists(ss_path):
    print(f"  sumstats already built: {os.path.basename(ss_path)}")
else:
    print(f"  Reading raw file ...")
    df = pd.read_csv(RAW_FILE, sep="\t", compression="gzip", low_memory=False)
    df = df.rename(columns={
        "rs_id":           "SNP",
        "effect_allele":   "A1",
        "other_allele":    "A2",
        "beta":            "BETA",
        "standard_error":  "SE",
    })
    ss = to_sumstats(df, N_SAMPLES)
    ss.to_csv(ss_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(ss):,} SNPs → {os.path.basename(ss_path)}")

python27 = get_python27()
print(f"\n── Running LDSC h2-cts (Path A) ──")
res_a = run_ldsc(ss_path, os.path.join(STUDY_DIR, f"{STUDY_ID}_ldsc_raw_direct"), python27)

# ── Path B: reuse existing harmonized result ──
print(f"\n── Path B: reusing {os.path.basename(RES_HM)} ──")
res_b = RES_HM

# ── Compare ──
if res_a and os.path.exists(res_b):
    r_a = pd.read_csv(res_a, sep="\t").sort_values("Coefficient_P_value")
    r_b = pd.read_csv(res_b, sep="\t").sort_values("Coefficient_P_value")
    _, r_a["FDR"] = fdrcorrection(r_a["Coefficient_P_value"].fillna(1))
    _, r_b["FDR"] = fdrcorrection(r_b["Coefficient_P_value"].fillna(1))

    merged = r_a[["Name","Coefficient_P_value"]].merge(
        r_b[["Name","Coefficient_P_value"]], on="Name", suffixes=("_a","_b")
    )
    rho, p = stats.spearmanr(merged["Coefficient_P_value_a"], merged["Coefficient_P_value_b"])

    print(f"\n{'='*60}")
    print(f"  RESULT — {STUDY_ID}  (inverse / Group B question)")
    print(f"{'='*60}")
    print(f"  Spearman ρ (222 cell types): {rho:.4f}  (p={p:.2e})")
    print(f"  Sig p<0.05  Path A (raw, no harmonisation):  {(r_a['Coefficient_P_value']<0.05).sum():>3}")
    print(f"  Sig p<0.05  Path B (raw → harmonizer.sh):    {(r_b['Coefficient_P_value']<0.05).sum():>3}")

    print(f"\n  Top 5 — Path A (raw → LDSC directly, no harmonisation):")
    print(r_a[["Name","Coefficient_P_value","FDR"]].head(5).to_string(index=False))
    print(f"\n  Top 5 — Path B (raw → harmonizer.sh → LDSC):")
    print(r_b[["Name","Coefficient_P_value","FDR"]].head(5).to_string(index=False))

    if rho > 0.95:
        verdict = "CONSISTENT — harmonisation does not change results on this raw file"
    elif rho > 0.80:
        verdict = "MOSTLY CONSISTENT — minor differences"
    else:
        verdict = "DISCORDANT — harmonisation changes results significantly on raw files"
    print(f"\n  Verdict: {verdict}")
else:
    print("ERROR: missing result file(s)")

print("\nDone.")
