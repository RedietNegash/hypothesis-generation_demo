#!/usr/bin/env python3
"""
run_multi_study_comparison.py
==============================
For each study: download OT pre-harmonised file (.h.tsv.gz), run it through
harmonizer.sh (Path A), run it through LDSC directly (Path B), compare ρ.

Studies (all hasSumstats=True in Open Targets):
  GCST90044364  Eczema/allergy (UKB)          N=455,449
  GCST006414    Atrial fibrillation            N=1,030,836
  GCST90474621  Height (UKB)                   N=457,377
  GCST90179150  BMI                            N=694,649
  GCST90444373  Alzheimer's disease            N=1,152,284
  GCST90018919  Schizophrenia                  N=629,347
"""
import os, json, subprocess, requests
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

# ── Constants ──────────────────────────────────────────────────────────────
OUT_DIR  = "/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation"
LDSC_WD  = "/mnt/hdd_1/rediet/hypothesis-generation-demo/notebooks"
W_HM3    = os.path.join(LDSC_WD, "data/OSF/w_hm3.snplist")
CTS_FILE = os.path.join(LDSC_WD, "data/PASS_AtrialFibrillation_Nielsen2018_sumstats_cell_types.cts")
BASELINE = "data/OSF/baseline_v1.2/baseline."
LDSCORES = "data/ldscores/all_merged_cCREs/all_merged_cCREs."
WEIGHTS  = "data/OSF/weights/weights.hm3_noMHC."
EBI_FTP  = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

HARMONIZER_SCRIPT  = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/harmonizer.sh"
HARMONIZER_REF_DIR = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/data/gwas_harm_ref"
HARMONIZER_REPO    = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser"
NEXTFLOW_ENV = {
    **os.environ,
    "PATH": "/mnt/hdd_1/rediet/hypothesis-generation-demo:/mnt/hdd_1/rediet/jdk-17/bin:" + os.environ.get("PATH", ""),
    "JAVA_HOME": "/mnt/hdd_1/rediet/jdk-17",
}
AMBIGUOUS = [{"A","T"}, {"C","G"}]

# ── Studies ────────────────────────────────────────────────────────────────
STUDIES = [
    {"id": "GCST90044364", "trait": "Eczema/allergy (UKB)",   "n": 455449},
    {"id": "GCST006414",   "trait": "Atrial fibrillation",    "n": 1030836},
    {"id": "GCST90474621", "trait": "Height (UKB)",           "n": 457377},
    {"id": "GCST90179150", "trait": "BMI",                    "n": 694649},
    {"id": "GCST90444373", "trait": "Alzheimer's disease",    "n": 1152284},
    {"id": "GCST90018919", "trait": "Schizophrenia",          "n": 629347},
]

# ── Helpers ────────────────────────────────────────────────────────────────

def gcst_range(study_id):
    num = int(study_id.replace("GCST", ""))
    lo  = (num // 1000) * 1000 + 1
    return f"GCST{lo:06d}-GCST{lo+999:06d}"


def find_hm_url(study_id):
    """Find the .h.tsv.gz URL on EBI FTP."""
    rng      = gcst_range(study_id)
    base_url = f"{EBI_FTP}/{rng}/{study_id}/harmonised"
    try:
        r = requests.get(base_url + "/", timeout=30)
        if r.status_code != 200:
            return None
        for line in r.text.split("\n"):
            if ".h.tsv.gz" in line and "build37" not in line.lower() and 'href="' in line:
                start = line.find('href="') + 6
                end   = line.find('"', start)
                fname = line[start:end]
                if fname.endswith(".h.tsv.gz"):
                    return f"{base_url}/{fname}"
    except Exception as e:
        print(f"  URL lookup error: {e}")
    return None


def download_file(url, dest):
    if os.path.exists(dest):
        print(f"  Already downloaded: {os.path.basename(dest)}")
        return True
    print(f"  Downloading: {os.path.basename(dest)} ...")
    try:
        r = requests.get(url, stream=True, timeout=120)
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)
        print(f"  Done: {os.path.getsize(dest)/1e6:.1f} MB")
        return True
    except Exception as e:
        print(f"  Download failed: {e}")
        return False


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


def reformat_for_harmonizer(hm_path, study_dir, study_id):
    out_dir  = os.path.join(study_dir, "gwas_reformatted_hm")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{study_id}.h.tsv.gz")
    if os.path.exists(out_path):
        print(f"  Already reformatted: {os.path.basename(out_path)}")
        return out_path
    df = pd.read_csv(hm_path, sep="\t", compression="gzip", low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cl = {c.lower(): c for c in df.columns}
    snp_col  = cl.get("hm_rsid") or cl.get("rsid") or cl.get("rs_id") or cl.get("snp")
    chr_col  = cl.get("hm_chrom") or cl.get("chromosome") or cl.get("chr")
    pos_col  = cl.get("hm_pos") or cl.get("base_pair_location") or cl.get("pos")
    a1_col   = cl.get("hm_effect_allele") or cl.get("effect_allele") or cl.get("a1")
    a2_col   = cl.get("hm_other_allele")  or cl.get("other_allele")  or cl.get("a2")
    beta_col = cl.get("hm_beta") or cl.get("beta") or cl.get("effect")
    se_col   = cl.get("standard_error") or cl.get("se")
    p_col    = cl.get("p_value") or cl.get("p-value") or cl.get("p")
    rename = {c: s for c, s in [(snp_col,"snp"),(chr_col,"chr"),(pos_col,"pos"),
                                  (a1_col,"a1"),(a2_col,"a2"),(beta_col,"beta"),
                                  (se_col,"se"),(p_col,"p")] if c}
    df = df.rename(columns=rename)
    if "chr" in df.columns:
        df["chr"] = df["chr"].astype(str).str.replace("chr","",regex=False)
    keep      = [c for c in ["chr","pos","snp","a1","a2","beta","se","p"] if c in df.columns]
    drop_cols = [c for c in ["a1","a2","beta","p"] if c in df.columns]
    df[keep].dropna(subset=drop_cols).to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(df):,} variants → {os.path.basename(out_path)}")
    return out_path


def run_harmonizer(reformatted_path, study_dir, study_id):
    ref_dir      = os.path.dirname(reformatted_path)
    expected_out = os.path.join(ref_dir, study_id)
    harm_dir     = os.path.join(study_dir, "harmonized_hm")
    os.makedirs(harm_dir, exist_ok=True)
    if os.path.exists(expected_out):
        print(f"  Harmonizer output already exists: {expected_out}")
        return expected_out
    log_file  = os.path.join(harm_dir, "harmonizer.log")
    input_abs = os.path.abspath(reformatted_path)
    orig_dir  = os.getcwd()
    os.chdir(harm_dir)
    try:
        print(f"  Running harmonizer.sh ...")
        with open(log_file, "w") as lf:
            subprocess.run([
                "bash", HARMONIZER_SCRIPT,
                "--input",     input_abs,
                "--build",     "GRCh38",
                "--ref",       HARMONIZER_REF_DIR,
                "--code-repo", HARMONIZER_REPO,
                "--threshold", "0.99",
            ], check=False, text=True, env=NEXTFLOW_ENV, stdout=lf, stderr=lf)
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
        ssf_path = reformatted_path if reformatted_path.endswith(".tsv.gz") \
                   else reformatted_path.replace(".gz", ".tsv.gz")
        print(f"  No final/ dir — SSF fallback: {os.path.basename(ssf_path)}")
        if os.path.exists(ssf_path):
            files     = [os.path.basename(ssf_path)]
            final_dir = os.path.dirname(ssf_path)
        else:
            print(f"  WARNING: SSF file not found: {ssf_path}")
            return None
    harmonized_file = os.path.join(final_dir, files[0])
    print(f"  Reading: {os.path.basename(harmonized_file)}")
    df = pd.read_csv(harmonized_file, sep="\t", compression="gzip", low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cl = {c.lower(): c for c in df.columns}
    snp_col  = cl.get("snp") or cl.get("hm_rsid") or cl.get("rsid")
    a1_col   = cl.get("a1") or cl.get("hm_effect_allele") or cl.get("effect_allele")
    a2_col   = cl.get("a2") or cl.get("hm_other_allele")  or cl.get("other_allele")
    beta_col = cl.get("beta") or cl.get("hm_beta")
    se_col   = cl.get("se") or cl.get("standard_error")
    rename = {c: s for c, s in [(snp_col,"SNP"),(a1_col,"A1"),(a2_col,"A2"),
                                  (beta_col,"BETA"),(se_col,"SE")] if c}
    df = df.rename(columns=rename)
    return to_sumstats(df, n)


def build_ot_sumstats(hm_path, n):
    """Path B: OT .h.tsv.gz → LDSC directly (no harmonizer)."""
    df = pd.read_csv(hm_path, sep="\t", compression="gzip", low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cl = {c.lower(): c for c in df.columns}
    snp_col  = cl.get("hm_rsid") or cl.get("rsid") or cl.get("rs_id") or cl.get("snp")
    a1_col   = cl.get("hm_effect_allele") or cl.get("effect_allele") or cl.get("a1")
    a2_col   = cl.get("hm_other_allele")  or cl.get("other_allele")  or cl.get("a2")
    beta_col = cl.get("hm_beta") or cl.get("beta")
    se_col   = cl.get("standard_error") or cl.get("se")
    rename = {c: s for c, s in [(snp_col,"SNP"),(a1_col,"A1"),(a2_col,"A2"),
                                  (beta_col,"BETA"),(se_col,"SE")] if c}
    df = df.rename(columns=rename)
    return to_sumstats(df, n)


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


def compare_results(res_a, res_b, study_id, trait):
    r_a = pd.read_csv(res_a, sep="\t").sort_values("Coefficient_P_value")
    r_b = pd.read_csv(res_b, sep="\t").sort_values("Coefficient_P_value")
    _, r_a["FDR"] = fdrcorrection(r_a["Coefficient_P_value"].fillna(1))
    _, r_b["FDR"] = fdrcorrection(r_b["Coefficient_P_value"].fillna(1))
    merged = r_a[["Name","Coefficient_P_value"]].merge(
        r_b[["Name","Coefficient_P_value"]], on="Name", suffixes=("_a","_b")
    )
    rho, p = stats.spearmanr(merged["Coefficient_P_value_a"], merged["Coefficient_P_value_b"])
    sig_a  = (r_a["Coefficient_P_value"] < 0.05).sum()
    sig_b  = (r_b["Coefficient_P_value"] < 0.05).sum()
    verdict = "CONSISTENT (ρ>0.95)" if rho > 0.95 else \
              "MOSTLY CONSISTENT"   if rho > 0.80 else "DISCORDANT"
    print(f"\n{'='*60}")
    print(f"  {study_id} — {trait}")
    print(f"  Spearman ρ: {rho:.4f}  (p={p:.2e})")
    print(f"  Sig p<0.05  Path A (via harmoniser): {sig_a:>3}")
    print(f"  Sig p<0.05  Path B (direct):         {sig_b:>3}")
    print(f"  Verdict: {verdict}")
    if sig_b > 0:
        print(f"\n  Top 5 significant cell types (Path B):")
        for _, row in r_b[r_b["Coefficient_P_value"] < 0.05].head(5).iterrows():
            print(f"    {row['Name']:<45} p={row['Coefficient_P_value']:.4f}")
    return {"study_id": study_id, "trait": trait, "rho": rho, "sig_a": sig_a, "sig_b": sig_b, "verdict": verdict}


# ── Main loop ──────────────────────────────────────────────────────────────
python27 = get_python27()
print(f"Python 2.7: {python27}\n")

summary = []
for s in STUDIES:
    sid   = s["id"]
    trait = s["trait"]
    n     = s["n"]
    study_dir = os.path.join(OUT_DIR, sid)
    os.makedirs(study_dir, exist_ok=True)

    print(f"\n{'#'*60}")
    print(f"  {sid} — {trait}  (N={n:,})")
    print(f"{'#'*60}")

    # Download OT harmonised file
    hm_file = os.path.join(study_dir, f"{sid}.h.tsv.gz")
    if not os.path.exists(hm_file):
        url = find_hm_url(sid)
        if not url:
            print(f"  ERROR: could not find harmonised URL for {sid} — skipping")
            continue
        print(f"  URL: {url}")
        if not download_file(url, hm_file):
            continue

    # ── Path A: .h.tsv.gz → harmonizer.sh → LDSC ──
    print("\n── Path A: .h.tsv.gz → harmonizer.sh → LDSC ──")
    res_a_path = os.path.join(study_dir, f"{sid}_ldsc_hm_via_harmonizer.cell_type_results.txt")
    if not os.path.exists(res_a_path):
        reformatted = reformat_for_harmonizer(hm_file, study_dir, sid)
        harm_out    = run_harmonizer(reformatted, study_dir, sid)
        ss_a        = harmonized_to_sumstats(harm_out, reformatted, n)
        if ss_a is None or len(ss_a) < 1000:
            print(f"  Too few SNPs — skipping {sid}")
            continue
        ss_a_path = os.path.join(study_dir, f"{sid}_hm_via_harmonizer.sumstats.gz")
        ss_a.to_csv(ss_a_path, sep="\t", index=False, compression="gzip")
        print(f"  Written {len(ss_a):,} SNPs → {os.path.basename(ss_a_path)}")
        res_a_path = run_ldsc(ss_a_path, os.path.join(study_dir, f"{sid}_ldsc_hm_via_harmonizer"), python27)
    else:
        print(f"  LDSC results already exist: {os.path.basename(res_a_path)}")

    # ── Path B: .h.tsv.gz → LDSC directly ──
    print("\n── Path B: .h.tsv.gz → LDSC directly ──")
    res_b_path = os.path.join(study_dir, f"{sid}_ldsc_ot.cell_type_results.txt")
    if not os.path.exists(res_b_path):
        ss_b = build_ot_sumstats(hm_file, n)
        if ss_b is None or len(ss_b) < 1000:
            print(f"  Too few SNPs — skipping {sid}")
            continue
        ss_b_path = os.path.join(study_dir, f"{sid}_ot.sumstats.gz")
        ss_b.to_csv(ss_b_path, sep="\t", index=False, compression="gzip")
        print(f"  Written {len(ss_b):,} SNPs → {os.path.basename(ss_b_path)}")
        res_b_path = run_ldsc(ss_b_path, os.path.join(study_dir, f"{sid}_ldsc_ot"), python27)
    else:
        print(f"  LDSC results already exist: {os.path.basename(res_b_path)}")

    # ── Compare ──
    if res_a_path and res_b_path and os.path.exists(res_a_path) and os.path.exists(res_b_path):
        result = compare_results(res_a_path, res_b_path, sid, trait)
        summary.append(result)

# ── Final summary ──────────────────────────────────────────────────────────
print(f"\n\n{'='*60}")
print("  FINAL SUMMARY")
print(f"{'='*60}")
print(f"  {'Study':<15} {'Trait':<30} {'ρ':>7}  {'Sig_B':>5}  Verdict")
print(f"  {'-'*70}")
for r in summary:
    print(f"  {r['study_id']:<15} {r['trait']:<30} {r['rho']:>7.4f}  {r['sig_b']:>5}  {r['verdict']}")

print("\nDone.")
