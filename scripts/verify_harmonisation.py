#!/usr/bin/env python3
import os, json, subprocess, requests
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

PARQUET      = "/mnt/hdd_1/rediet/ad-regulome/data/opentargets_studies/studies.parquet"
OUT_DIR      = "/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation"
LDSC_WD      = "/mnt/hdd_1/rediet/hypothesis-generation-demo/notebooks"
W_HM3        = "data/OSF/w_hm3.snplist"
CTS_FILE     = "data/PASS_AtrialFibrillation_Nielsen2018_sumstats_cell_types.cts"
BASELINE     = "data/OSF/baseline_v1.2/baseline."
LDSCORES     = "data/ldscores/all_merged_cCREs/all_merged_cCREs."
WEIGHTS      = "data/OSF/weights/weights.hm3_noMHC."
EBI_FTP      = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
N_STUDIES    = 20
RANDOM_SEED  = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Step 1: Sampling studies from Open Targets parquet")
print("=" * 60)

df = pd.read_parquet(PARQUET)

eligible = df[
    (df["studyType"] == "gwas") &
    (df["hasSumstats"].fillna(False) == True) &
    (
        (df["projectId"] == "FINNGEN_R12") |
        df["sumstatQCValues"].apply(
            lambda x: isinstance(x, (list, np.ndarray)) and len(x) > 0
        )
    )
][["studyId", "projectId", "traitFromSource", "nSamples", "nCases",
   "nControls", "publicationFirstAuthor", "publicationDate"]].copy()

sample = eligible.sample(N_STUDIES, random_state=RANDOM_SEED).reset_index(drop=True)
sample["trait"] = sample["traitFromSource"].str[:50]

print(f"  Total eligible studies : {len(eligible):,}")
print(f"    FinnGen R12          : {(eligible['projectId']=='FINNGEN_R12').sum():,}")
print(f"    GCST with QC values  : {(eligible['projectId']!='FINNGEN_R12').sum():,}")
print(f"\n  Random sample of {N_STUDIES}:\n")
print(sample[["studyId", "projectId", "trait", "nSamples", "nCases"]].to_string(index=False))


def gcst_range(study_id):
    num = int(study_id.replace("GCST", ""))
    lo  = (num // 1000) * 1000 + 1
    hi  = lo + 999
    return f"GCST{lo:06d}-GCST{hi:06d}"


def find_harmonised_url(study_id):
    rng  = gcst_range(study_id)
    url  = f"{EBI_FTP}/{rng}/{study_id}/harmonised/"
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        return None
    for line in resp.text.split("\n"):
        if ".h.tsv.gz" in line and "build37" not in line:
            start = line.find('href="') + 6
            end   = line.find('"', start)
            fname = line[start:end]
            if fname.endswith(".h.tsv.gz"):
                return url + fname
    return None


def download_file(url, out_path):
    if os.path.exists(out_path):
        print(f"    Already downloaded: {os.path.basename(out_path)}")
        return True
    print(f"    Downloading {url} ...")
    try:
        resp = requests.get(url, stream=True, timeout=120)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=4 * 1024 * 1024):
                f.write(chunk)
        print(f"    Saved: {os.path.basename(out_path)}")
        return True
    except Exception as e:
        print(f"    ERROR downloading: {e}")
        return False


AMBIGUOUS = [{"A", "T"}, {"C", "G"}]

def load_hm3():
    path = os.path.join(LDSC_WD, W_HM3)
    return pd.read_csv(path, sep="\t")[["SNP", "A1", "A2"]]

def remove_ambiguous(df):
    return df[~df.apply(lambda r: set([r["A1"], r["A2"]]) in AMBIGUOUS, axis=1)]


def build_sumstats_harmonised(raw_path, n_samples):
    needed = ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error"]
    df = pd.read_csv(raw_path, sep="\t", compression="gzip", usecols=needed, low_memory=False)
    df = df.rename(columns={
        "hm_rsid":          "SNP",
        "hm_effect_allele": "A1",
        "hm_other_allele":  "A2",
        "hm_beta":          "BETA",
        "standard_error":   "SE",
    })
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"]   = pd.to_numeric(df["SE"],   errors="coerce")
    df["A1"]   = df["A1"].str.upper()
    df["A2"]   = df["A2"].str.upper()
    df["Z"]    = df["BETA"] / df["SE"]
    df["N"]    = n_samples
    df = df.dropna(subset=["SNP", "Z"])
    hm3 = load_hm3()
    df  = df.merge(hm3, on="SNP", suffixes=("", "_hm3"))
    df  = df[(df["A1"] == df["A1_hm3"]) | (df["A1"] == df["A2_hm3"])]
    df  = df.drop(columns=["A1_hm3", "A2_hm3"])
    df  = remove_ambiguous(df)
    return df[["SNP", "A1", "A2", "Z", "N"]].dropna()


def build_sumstats_raw(raw_path, n_samples):
    peek = pd.read_csv(raw_path, sep="\t", compression="gzip", nrows=0, low_memory=False)
    cols = peek.columns.tolist()

    snp_col  = next((c for c in ["hm_rsid", "rsid", "variant_id", "snp", "SNP"] if c in cols), None)
    a1_col   = next((c for c in ["effect_allele", "a1", "A1", "alt"]            if c in cols), None)
    a2_col   = next((c for c in ["other_allele",  "a2", "A2", "ref"]            if c in cols), None)
    beta_col = next((c for c in ["beta", "BETA", "effect", "log_odds"]          if c in cols), None)
    se_col   = next((c for c in ["standard_error", "se", "SE", "sebeta"]        if c in cols), None)

    if not all([snp_col, a1_col, a2_col, beta_col, se_col]):
        print(f"    WARNING: could not find raw columns. Available: {cols[:10]}")
        return None

    df = pd.read_csv(raw_path, sep="\t", compression="gzip",
                     usecols=[snp_col, a1_col, a2_col, beta_col, se_col], low_memory=False)
    df = df.rename(columns={snp_col: "SNP", a1_col: "A1", a2_col: "A2",
                             beta_col: "BETA", se_col: "SE"})
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"]   = pd.to_numeric(df["SE"],   errors="coerce")
    df["A1"]   = df["A1"].str.upper()
    df["A2"]   = df["A2"].str.upper()
    df["Z"]    = df["BETA"] / df["SE"]
    df["N"]    = n_samples
    df = df.dropna(subset=["SNP", "Z"])
    hm3 = load_hm3()
    df  = df.merge(hm3, on="SNP", suffixes=("", "_hm3"))
    df  = df[(df["A1"] == df["A1_hm3"]) | (df["A1"] == df["A2_hm3"])]
    df  = df.drop(columns=["A1_hm3", "A2_hm3"])
    df  = remove_ambiguous(df)
    return df[["SNP", "A1", "A2", "Z", "N"]].dropna()


def get_python27():
    j    = subprocess.run(["conda", "env", "list", "--json"], capture_output=True, text=True, check=True)
    envs = json.loads(j.stdout)["envs"]
    return os.path.join([e for e in envs if "ldsc27" in e][0], "bin", "python")


def run_ldsc(sumstats_path, out_prefix, python27):
    results_file = f"{out_prefix}.cell_type_results.txt"
    if os.path.exists(results_file):
        print(f"    Results already exist: {os.path.basename(results_file)}")
        return results_file
    cmd = [
        python27, os.path.join(LDSC_WD, "tools/ldsc/ldsc.py"),
        "--h2-cts",         sumstats_path,
        "--ref-ld-chr",     f"{BASELINE},{LDSCORES}",
        "--ref-ld-chr-cts", os.path.join(LDSC_WD, CTS_FILE),
        "--w-ld-chr",       os.path.join(LDSC_WD, WEIGHTS),
        "--out",            out_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=LDSC_WD)
    if result.returncode != 0:
        print(f"    LDSC ERROR:\n{result.stderr[-500:]}")
        return None
    return results_file


def compare_results(study_id, res_hm, res_raw):
    r_hm  = pd.read_csv(res_hm,  sep="\t").sort_values("Coefficient_P_value")
    r_raw = pd.read_csv(res_raw, sep="\t").sort_values("Coefficient_P_value")
    _, r_hm["FDR"]  = fdrcorrection(r_hm["Coefficient_P_value"].fillna(1))
    _, r_raw["FDR"] = fdrcorrection(r_raw["Coefficient_P_value"].fillna(1))

    merged = r_hm[["Name", "Coefficient_P_value"]].merge(
        r_raw[["Name", "Coefficient_P_value"]], on="Name", suffixes=("_hm", "_raw")
    )
    rho, p = stats.spearmanr(merged["Coefficient_P_value_hm"], merged["Coefficient_P_value_raw"])

    print(f"\n  {'─'*54}")
    print(f"  {study_id}")
    print(f"  {'─'*54}")
    print(f"  Rank correlation (Spearman ρ): {rho:.4f}  (p={p:.2e})")
    print(f"  Sig p<0.05  — HM: {(r_hm['Coefficient_P_value']<0.05).sum():>3}  |  Raw: {(r_raw['Coefficient_P_value']<0.05).sum():>3}")
    print(f"  Sig FDR<0.05— HM: {(r_hm['FDR']<0.05).sum():>3}  |  Raw: {(r_raw['FDR']<0.05).sum():>3}")
    print(f"\n  Top 5 — WITH harmonisation:")
    print(r_hm[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))
    print(f"\n  Top 5 — WITHOUT harmonisation (raw):")
    print(r_raw[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))

    if rho > 0.95:
        verdict = "✅ CONSISTENT — harmonisation step not needed"
    elif rho > 0.80:
        verdict = "⚠️  MOSTLY CONSISTENT — minor differences"
    else:
        verdict = "❌ DISCORDANT — harmonisation changes results significantly"
    print(f"\n  Verdict: {verdict}")
    return {"study_id": study_id, "spearman_rho": rho, "p": p,
            "sig_hm": (r_hm["FDR"]<0.05).sum(), "sig_raw": (r_raw["FDR"]<0.05).sum(),
            "verdict": verdict}


python27 = get_python27()
print(f"\nPython 2.7 for LDSC: {python27}\n")

summary = []
skipped = []

for _, row in sample.iterrows():
    study_id = row["studyId"]
    project  = row["projectId"]
    n        = int(row["nSamples"]) if pd.notna(row["nSamples"]) else 100000

    print(f"\n{'#'*60}")
    print(f"  {study_id}  |  {row['trait']}")
    print(f"  project={project}  N={n:,}  cases={row['nCases']}")
    print(f"{'#'*60}")

    if project == "FINNGEN_R12":
        print("  FinnGen file is on GCS — requires gsutil + credentials, skipping.")
        skipped.append({"study_id": study_id, "reason": "FinnGen GCS"})
        continue

    print("  Locating harmonised file on EBI FTP ...")
    url = find_harmonised_url(study_id)
    if not url:
        print("  No harmonised file found — skipping.")
        skipped.append({"study_id": study_id, "reason": "No harmonised file on EBI FTP"})
        continue
    print(f"  URL: {url}")

    raw_path = os.path.join(OUT_DIR, f"{study_id}.h.tsv.gz")
    if not download_file(url, raw_path):
        skipped.append({"study_id": study_id, "reason": "Download failed"})
        continue

    print("  Building sumstats WITH harmonisation (hm_* columns) ...")
    ss_hm = build_sumstats_harmonised(raw_path, n)
    if ss_hm is None or len(ss_hm) < 1000:
        print(f"  Too few SNPs ({len(ss_hm) if ss_hm is not None else 0}) — skipping.")
        skipped.append({"study_id": study_id, "reason": "Too few SNPs (harmonised)"})
        continue
    ss_hm_path = os.path.join(OUT_DIR, f"{study_id}_harmonised.sumstats.gz")
    ss_hm.to_csv(ss_hm_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(ss_hm):,} SNPs → {os.path.basename(ss_hm_path)}")

    print("  Building sumstats WITHOUT harmonisation (raw columns) ...")
    ss_raw = build_sumstats_raw(raw_path, n)
    if ss_raw is None or len(ss_raw) < 1000:
        print(f"  Too few SNPs ({len(ss_raw) if ss_raw is not None else 0}) — skipping.")
        skipped.append({"study_id": study_id, "reason": "Too few SNPs (raw)"})
        continue
    ss_raw_path = os.path.join(OUT_DIR, f"{study_id}_raw.sumstats.gz")
    ss_raw.to_csv(ss_raw_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(ss_raw):,} SNPs → {os.path.basename(ss_raw_path)}")

    shared = set(ss_hm["SNP"]) & set(ss_raw["SNP"])
    z_hm   = ss_hm.set_index("SNP").loc[list(shared), "Z"]
    z_raw  = ss_raw.set_index("SNP").loc[list(shared), "Z"]
    r, _   = stats.pearsonr(z_hm, z_raw)
    print(f"  Z-score Pearson r (shared SNPs, n={len(shared):,}): {r:.4f}")

    print("  Running LDSC h2-cts (harmonised) ...")
    res_hm  = run_ldsc(ss_hm_path,  os.path.join(OUT_DIR, f"{study_id}_ldsc_harmonised"), python27)

    print("  Running LDSC h2-cts (raw) ...")
    res_raw = run_ldsc(ss_raw_path, os.path.join(OUT_DIR, f"{study_id}_ldsc_raw"), python27)

    if res_hm and res_raw:
        result = compare_results(study_id, res_hm, res_raw)
        result["z_pearson_r"] = r
        summary.append(result)
    else:
        skipped.append({"study_id": study_id, "reason": "LDSC failed"})


print(f"\n\n{'='*60}")
print("  FINAL SUMMARY")
print(f"{'='*60}")

if summary:
    s = pd.DataFrame(summary)
    print(s[["study_id", "z_pearson_r", "spearman_rho", "sig_hm", "sig_raw", "verdict"]].to_string(index=False))
    print(f"\n  Studies analysed     : {len(s)}")
    print(f"  Consistent (ρ>0.95)  : {(s['spearman_rho'] > 0.95).sum()} / {len(s)}")

if skipped:
    print(f"\n  Skipped ({len(skipped)}):")
    for sk in skipped:
        print(f"    {sk['study_id']}: {sk['reason']}")

print("\nDone.")
