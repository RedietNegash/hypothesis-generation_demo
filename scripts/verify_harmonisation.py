#!/usr/bin/env python3
"""
verify_harmonisation.py
=======================
Tests two hypotheses with LDSC h2-cts (222 cell types):

  Group A (pre-harmonised) — FinnGen R12 + GCST with sumstatQCValues
    Path A: OT/FinnGen file directly → LDSC   (skip harmoniser)
    Path B: Raw file → harmonizer.sh → LDSC
    Expect: ρ > 0.95  →  safe to skip harmonisation

  Group B (inverse / not harmonised) — GCST without sumstatQCValues
    Path A: Raw file → LDSC                    (no harmonisation)
    Path B: Raw file → harmonizer.sh → LDSC
    Expect: ρ < 0.95  →  harmonisation matters for un-harmonised studies
"""
import os, json, glob, subprocess, requests
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

pd.set_option("future.no_silent_downcasting", True)

PARQUET      = "/mnt/hdd_1/rediet/ad-regulome/data/opentargets_studies/studies.parquet"
OUT_DIR      = "/mnt/hdd_1/rediet/ad-regulome/data/gwas/verify_harmonisation"
LDSC_WD      = "/mnt/hdd_1/rediet/hypothesis-generation-demo/notebooks"
W_HM3        = os.path.join(LDSC_WD, "data/OSF/w_hm3.snplist")
CTS_FILE     = os.path.join(LDSC_WD, "data/PASS_AtrialFibrillation_Nielsen2018_sumstats_cell_types.cts")
BASELINE     = "data/OSF/baseline_v1.2/baseline."
LDSCORES     = "data/ldscores/all_merged_cCREs/all_merged_cCREs."
WEIGHTS      = "data/OSF/weights/weights.hm3_noMHC."
EBI_FTP      = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
HM_LIST_URL  = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt"

HARMONIZER_SCRIPT   = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/harmonizer.sh"
HARMONIZER_REF_DIR  = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/data/gwas_harm_ref"
HARMONIZER_REPO     = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser"
NEXTFLOW_ENV = {
    **os.environ,
    "PATH": "/mnt/hdd_1/rediet/hypothesis-generation-demo:/mnt/hdd_1/rediet/jdk-17/bin:" + os.environ.get("PATH", ""),
    "JAVA_HOME": "/mnt/hdd_1/rediet/jdk-17",
}

N_STUDIES    = 1      # studies to sample per group (A and B)
RANDOM_SEED  = 42
AMBIGUOUS    = [{"A", "T"}, {"C", "G"}]

os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: Build harmonised-URL lookup (one download replaces per-study FTP scraping) ──
print("=" * 60)
print("Step 1: Building harmonised URL lookup from EBI list")
print("=" * 60)


def build_hm_lookup():
    """
    Download harmonised_list.txt once (cached) → {study_id: full_harmonised_url}.
    Replaces the old approach of making 2 HTTP requests per study to scrape
    FTP directory listings.  Format of each line:
      ./GCST90615001-GCST90616000/GCST90615056/harmonised/GCST90615056.h.tsv.gz
    """
    cache = os.path.join(OUT_DIR, "harmonised_list.txt")
    if not os.path.exists(cache):
        print(f"  Downloading {HM_LIST_URL} ...")
        resp = requests.get(HM_LIST_URL, timeout=300)
        resp.raise_for_status()
        with open(cache, "w") as f:
            f.write(resp.text)
        print(f"  Saved → {cache}")
    else:
        print(f"  Using cached harmonised_list.txt")

    lookup = {}
    with open(cache) as f:
        for line in f:
            line = line.strip()
            if not line.endswith(".h.tsv.gz"):
                continue
            parts = line.lstrip("./").split("/")
            if len(parts) < 4:
                continue
            rng_folder, study_id, _, fname = parts[0], parts[1], parts[2], parts[-1]
            if "build37" in fname.lower():
                continue
            lookup[study_id] = f"{EBI_FTP}/{rng_folder}/{study_id}/harmonised/{fname}"

    print(f"  {len(lookup):,} harmonised studies indexed")
    return lookup


HM_LOOKUP = build_hm_lookup()

# ── Step 2: Sample studies — Group A (harmonised) + Group B (inverse) ──────
print("\n" + "=" * 60)
print("Step 2: Sampling studies from Open Targets parquet")
print("=" * 60)

df = pd.read_parquet(PARQUET)

has_qc = df["sumstatQCValues"].apply(
    lambda x: isinstance(x, (list, np.ndarray)) and len(x) > 0
)
base_filter = (df["studyType"] == "gwas") & df["hasSumstats"].fillna(False)
_cols = ["studyId", "projectId", "traitFromSource", "nSamples", "nCases",
         "nControls", "publicationFirstAuthor", "summarystatsLocation"]

# Group A: FinnGen R12 OR GCST with QC values (OT already harmonised these)
eligible_a = df[base_filter & ((df["projectId"] == "FINNGEN_R12") | has_qc)][_cols].copy()
eligible_a["group"] = "A"

# Group B (inverse): GCST studies WITHOUT sumstatQCValues — not harmonised by OT
eligible_b = df[base_filter & (df["projectId"] == "GCST") & ~has_qc][_cols].copy()
eligible_b["group"] = "B"

sample_a = eligible_a.sample(N_STUDIES, random_state=RANDOM_SEED).reset_index(drop=True)
sample_b = eligible_b.sample(N_STUDIES, random_state=RANDOM_SEED).reset_index(drop=True)
sample   = pd.concat([sample_a, sample_b], ignore_index=True)
sample["trait"] = sample["traitFromSource"].str[:45]

print(f"  Group A eligible : {len(eligible_a):,}  "
      f"(FinnGen={( eligible_a['projectId']=='FINNGEN_R12').sum():,}  "
      f"GCST={( eligible_a['projectId']!='FINNGEN_R12').sum():,})")
print(f"  Group B eligible : {len(eligible_b):,}  (GCST without QC values)")
print(f"\n  Sample ({N_STUDIES} per group):\n")
print(sample[["studyId", "projectId", "group", "trait", "nSamples"]].to_string(index=False))


def gcst_range(study_id):
    num = int(study_id.replace("GCST", ""))
    lo  = (num // 1000) * 1000 + 1
    return f"GCST{lo:06d}-GCST{lo+999:06d}"


def find_raw_url(study_id):
    """
    Get the raw summary-stats URL from the EBI FTP root directory listing.
    One HTTP request per study; harmonised URL now comes from HM_LOOKUP.
    """
    base_url = f"{EBI_FTP}/{gcst_range(study_id)}/{study_id}"
    try:
        resp = requests.get(base_url + "/", timeout=30)
        if resp.status_code == 200:
            for line in resp.text.split("\n"):
                if 'href="' in line and "Parent" not in line:
                    s = line.find('href="') + 6
                    e = line.find('"', s)
                    fname = line[s:e]
                    if fname.endswith((".gz", ".txt", ".tsv")) and "/" not in fname:
                        return f"{base_url}/{fname}"
    except Exception as ex:
        print(f"    WARNING: FTP listing failed: {ex}")
    return None


def download_finngen(gcs_path, out_path):
    """Download a FinnGen file from GCS using the summarystatsLocation column."""
    if os.path.exists(out_path):
        print(f"    Already downloaded: {os.path.basename(out_path)}")
        return True
    print(f"    gsutil cp {gcs_path} {out_path}")
    result = subprocess.run(
        ["gsutil", "cp", gcs_path, out_path],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode == 0:
        print(f"    Saved: {os.path.basename(out_path)}")
        return True
    print(f"    gsutil ERROR:\n{result.stderr[-300:]}")
    return False


def download_file(url, out_path):
    if os.path.exists(out_path):
        print(f"    Already downloaded: {os.path.basename(out_path)}")
        return True
    print(f"    Downloading {url} ...")
    try:
        resp = requests.get(url, stream=True, timeout=300)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=4 * 1024 * 1024):
                f.write(chunk)
        print(f"    Saved: {os.path.basename(out_path)}")
        return True
    except Exception as e:
        print(f"    ERROR downloading: {e}")
        return False


def load_hm3():
    return pd.read_csv(W_HM3, sep="\t")[["SNP", "A1", "A2"]]


def to_sumstats(df, n_samples):
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
    df  = df[~df.apply(lambda r: set([r["A1"], r["A2"]]) in AMBIGUOUS, axis=1)]
    return df[["SNP", "A1", "A2", "Z", "N"]].dropna()


def build_sumstats_ot(hm_path, n_samples):
    """WITHOUT harmonisation: read OT harmonised file.
    Handles two EBI formats:
      - hm_ prefix columns (hm_rsid, hm_beta, hm_effect_allele, hm_other_allele)
      - plain columns (rsid/rs_id, beta, effect_allele, other_allele)
    """
    peek = pd.read_csv(hm_path, sep="\t", compression="gzip", nrows=0, low_memory=False)
    cols = set(peek.columns.str.lower())
    print(f"    Harmonised file columns: {list(peek.columns)}")

    if "hm_rsid" in cols:
        needed = ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error"]
        df = pd.read_csv(hm_path, sep="\t", compression="gzip", usecols=needed, low_memory=False)
        df = df.rename(columns={
            "hm_rsid":          "SNP",
            "hm_effect_allele": "A1",
            "hm_other_allele":  "A2",
            "hm_beta":          "BETA",
            "standard_error":   "SE",
        })
    else:
        cl = {c.lower(): c for c in peek.columns}
        snp_col  = cl.get("hm_rsid") or cl.get("rsid") or cl.get("rs_id") or cl.get("snp") or cl.get("variant_id")
        a1_col   = cl.get("hm_effect_allele") or cl.get("effect_allele") or cl.get("a1")
        a2_col   = cl.get("hm_other_allele")  or cl.get("other_allele")  or cl.get("a2")
        beta_col = cl.get("hm_beta") or cl.get("beta") or cl.get("effect")
        se_col   = cl.get("standard_error") or cl.get("se") or cl.get("stderr")

        if not all([snp_col, a1_col, a2_col, beta_col, se_col]):
            print(f"    WARNING: could not map harmonised columns. Available: {list(peek.columns)}")
            return None

        print(f"    Mapped: SNP={snp_col} A1={a1_col} A2={a2_col} BETA={beta_col} SE={se_col}")
        needed = [snp_col, a1_col, a2_col, beta_col, se_col]
        df = pd.read_csv(hm_path, sep="\t", compression="gzip", usecols=needed, low_memory=False)
        df = df.rename(columns={snp_col: "SNP", a1_col: "A1", a2_col: "A2",
                                 beta_col: "BETA", se_col: "SE"})

    return to_sumstats(df, n_samples)


def build_sumstats_raw(reformatted_path, n_samples):
    """
    Group B — Path A: use the reformatted raw file directly (no harmonisation).
    reformat_gwas() writes lowercase SSF column names (snp, a1, a2, beta, se),
    so we just rename them to the uppercase format to_sumstats() expects.
    """
    df = pd.read_csv(reformatted_path, sep="\t", compression="gzip", low_memory=False)
    df = df.rename(columns={"snp": "SNP", "a1": "A1", "a2": "A2",
                             "beta": "BETA", "se": "SE"})
    required = {"SNP", "A1", "A2", "BETA", "SE"}
    if not required.issubset(df.columns):
        print(f"    Missing columns {required - set(df.columns)}. Have: {list(df.columns)}")
        return None
    print(f"    Read {len(df):,} variants from raw reformatted file (no harmonisation)")
    return to_sumstats(df, n_samples)


def reformat_gwas(raw_path, study_dir):
    """
    Step 3 from the AD regulome notebook:
    Auto-detect columns in raw file, rename to standard names,
    write to <study_dir>/gwas_reformatted/<file>.gz
    """
    import gzip as _gzip

    out_dir = os.path.join(study_dir, "gwas_reformatted")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, os.path.basename(raw_path))

    if os.path.exists(out_path):
        print(f"    Reformatted file already exists: {os.path.basename(out_path)}")
        return out_path

    _lower   = raw_path.lower()
    _is_gz   = _lower.endswith(".gz") or _lower.endswith(".bgz")
    _open_fn = _gzip.open if _is_gz else open
    _comp    = "gzip" if _is_gz else None

    with _open_fn(raw_path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                _sep = "\t" if "\t" in line else " "
                break

    df = pd.read_csv(raw_path, sep=_sep, compression=_comp, low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    print(f"    Raw columns: {list(df.columns)}")

    cl = {c.lower(): c for c in df.columns}

    snp_col  = cl.get("markername") or cl.get("snptestid") or cl.get("id") or cl.get("snp") or cl.get("rsid") or cl.get("variant_id") or cl.get("rs_id")
    chr_col  = cl.get("chr") or cl.get("chrom") or cl.get("chromosome") or cl.get("#chrom") or cl.get("hm_chrom")
    pos_col  = cl.get("pos") or cl.get("bp") or cl.get("position") or cl.get("bp_hg19") or cl.get("base_pair_location")
    a1_col   = cl.get("a1") or cl.get("effect_allele") or cl.get("alt") or cl.get("allele1") or cl.get("hm_effect_allele")
    a2_col   = cl.get("a2") or cl.get("noneffect_allele") or cl.get("other_allele") or cl.get("ref") or cl.get("allele2") or cl.get("hm_other_allele")
    beta_col = cl.get("logor") or cl.get("log_or") or cl.get("beta") or cl.get("b") or cl.get("effect") or cl.get("hm_beta")
    se_col   = cl.get("se_gc") or cl.get("stderrlogor") or cl.get("se") or cl.get("stderr") or cl.get("standard_error") or cl.get("sebeta")
    p_col    = cl.get("p-value_gc") or cl.get("pvalue") or cl.get("p_value") or cl.get("p-value") or cl.get("p")
    n_col    = cl.get("n_samples") or cl.get("neff") or cl.get("n") or cl.get("n_total")

    rename = {}
    if snp_col:  rename[snp_col]  = "snp"
    if a1_col:   rename[a1_col]   = "a1"
    if a2_col:   rename[a2_col]   = "a2"
    if beta_col: rename[beta_col] = "beta"
    if se_col:   rename[se_col]   = "se"
    if p_col:    rename[p_col]    = "p"
    if n_col:    rename[n_col]    = "n"
    if chr_col:  rename[chr_col]  = "chr"
    if pos_col:  rename[pos_col]  = "pos"
    df = df.rename(columns=rename)

    if "chr" not in df.columns or "pos" not in df.columns:
        print("    No chr/pos — looking up from bim files ...")
        bim_files = sorted(glob.glob(
            os.path.join(LDSC_WD, "data/reference/GRCh38/plink_files/1000G.EUR.hg38.*.bim")
        ))
        if bim_files:
            bim_map = pd.concat([
                pd.read_csv(f, sep="\t", header=None, names=["chr","snp","cm","pos","a1b","a2b"])[["snp","chr","pos"]]
                for f in bim_files
            ]).drop_duplicates("snp").set_index("snp")
            df["chr"] = df["snp"].map(bim_map["chr"])
            df["pos"] = df["snp"].map(bim_map["pos"])
            df = df.dropna(subset=["chr","pos"])
            df["chr"] = df["chr"].astype(int).astype(str)
            df["pos"] = df["pos"].astype(int).astype(str)
            print(f"    Mapped {len(df):,} variants with chr/pos")

    if "chr" in df.columns:
        df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)

    keep = [c for c in ["chr","pos","snp","a1","a2","beta","se","p","n"] if c in df.columns]
    dropna_cols = [c for c in ["a1","a2","beta","p"] if c in df.columns]
    df[keep].dropna(subset=dropna_cols).to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"    Written {len(df):,} variants to {os.path.basename(out_path)}")
    return out_path


def run_harmonizer(reformatted_path, study_dir):
    """
    Step 5 from the AD regulome notebook:
    Run harmonizer.sh on the reformatted file.
    Returns path to the harmonized output directory, or None on failure.

    The harmonizer creates output in <gwas_reformatted_dir>/<study_id>/
    If no final/ is present (file already GRCh38), the SSF file itself is the output.
    """
    study_id     = os.path.basename(study_dir)
    ref_dir      = os.path.dirname(reformatted_path)   # .../gwas_reformatted/
    expected_out = os.path.join(ref_dir, study_id)     # .../gwas_reformatted/<study_id>/
    harm_work_dir = os.path.join(study_dir, "harmonized")
    os.makedirs(harm_work_dir, exist_ok=True)

    if os.path.exists(expected_out):
        print(f"    Harmonizer output already exists: {expected_out}")
        return expected_out

    log_file  = os.path.join(harm_work_dir, "harmonizer.log")
    input_abs = os.path.abspath(reformatted_path)
    orig_dir  = os.getcwd()

    os.chdir(harm_work_dir)
    try:
        print(f"    Running harmonizer.sh on {os.path.basename(reformatted_path)} ...")
        with open(log_file, "w") as lf:
            result = subprocess.run([
                "bash", HARMONIZER_SCRIPT,
                "--input",     input_abs,
                "--build",     "GRCh38",
                "--ref",       HARMONIZER_REF_DIR,
                "--code-repo", HARMONIZER_REPO,
                "--threshold", "0.99",
            ], check=False, text=True, env=NEXTFLOW_ENV, stdout=lf, stderr=lf)

        print(f"    harmonizer.sh exit code: {result.returncode}")

        if os.path.exists(expected_out):
            print(f"    Harmonized output dir: {expected_out}")
            return expected_out
        else:
            print(f"    No study subdir found — returning ref_dir for SSF fallback")
            return ref_dir

    except Exception as e:
        print(f"    ERROR running harmonizer: {e}")
        return None
    finally:
        os.chdir(orig_dir)


def harmonized_to_sumstats(harmonized_output_dir, reformatted_path, n_samples):
    """
    Step 6 from the AD regulome notebook:
    Read harmonizer output, rename columns, compute Z, merge HM3, remove ambiguous.
    Falls back to SSF file if no final/ directory (file already in GRCh38).
    Returns a sumstats DataFrame.
    """
    final_dir = os.path.join(harmonized_output_dir, "final")
    files = []
    if os.path.exists(final_dir):
        files = [f for f in os.listdir(final_dir) if f.endswith(".tsv.gz")]

    if not files:
        # Notebook fallback: use SSF file directly when no final/ dir exists
        ssf_path = reformatted_path.replace(".gz", ".tsv.gz")
        print(f"    No final/ dir — trying SSF fallback: {os.path.basename(ssf_path)}")
        if os.path.exists(ssf_path):
            files     = [os.path.basename(ssf_path)]
            final_dir = os.path.dirname(ssf_path)
        else:
            print(f"    WARNING: SSF file not found: {ssf_path}")
            return None

    harmonized_file = os.path.join(final_dir, files[0])
    print(f"    Reading harmonized output: {os.path.basename(harmonized_file)}")
    df = pd.read_csv(harmonized_file, sep="\t", compression="gzip")
    print(f"    Columns: {list(df.columns)}")

    df = df.rename(columns={
        "rsid":           "SNP",
        "effect_allele":  "A1",
        "other_allele":   "A2",
        "beta":           "BETA",
        "standard_error": "SE",
        "p_value":        "P",
    })

    if "N" not in df.columns:
        df["N"] = n_samples

    return to_sumstats(df, n_samples)


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
        "--ref-ld-chr-cts", CTS_FILE,
        "--w-ld-chr",       os.path.join(LDSC_WD, WEIGHTS),
        "--out",            out_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=LDSC_WD)
    if result.returncode != 0:
        print(f"    LDSC ERROR:\n{result.stderr[-500:]}")
        return None
    return results_file


def compare_results(study_id, res_a, res_b, group):
    """Compare LDSC h2-cts results from Path A vs Path B for either group."""
    r_a = pd.read_csv(res_a, sep="\t").sort_values("Coefficient_P_value")
    r_b = pd.read_csv(res_b, sep="\t").sort_values("Coefficient_P_value")
    _, r_a["FDR"] = fdrcorrection(r_a["Coefficient_P_value"].fillna(1))
    _, r_b["FDR"] = fdrcorrection(r_b["Coefficient_P_value"].fillna(1))

    merged = r_a[["Name", "Coefficient_P_value"]].merge(
        r_b[["Name", "Coefficient_P_value"]], on="Name", suffixes=("_a", "_b")
    )
    rho, p = stats.spearmanr(merged["Coefficient_P_value_a"], merged["Coefficient_P_value_b"])

    if group == "A":
        label_a = "OT harmonised  (skip harmoniser)"
        label_b = "Raw → harmonizer.sh"
    else:
        label_a = "Raw file       (no harmonisation)"
        label_b = "Raw → harmonizer.sh"

    print(f"\n  {'─'*54}")
    print(f"  {study_id}  [Group {group}]")
    print(f"  {'─'*54}")
    print(f"  Spearman ρ (222 cell types): {rho:.4f}  (p={p:.2e})")
    print(f"  Sig p<0.05   Path A: {(r_a['Coefficient_P_value']<0.05).sum():>3}  |  Path B: {(r_b['Coefficient_P_value']<0.05).sum():>3}")
    print(f"  Sig FDR<0.05 Path A: {(r_a['FDR']<0.05).sum():>3}  |  Path B: {(r_b['FDR']<0.05).sum():>3}")
    print(f"\n  Top 5 — Path A ({label_a}):")
    print(r_a[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))
    print(f"\n  Top 5 — Path B ({label_b}):")
    print(r_b[["Name", "Coefficient_P_value", "FDR"]].head(5).to_string(index=False))

    if group == "A":
        # Harmonised studies: expect high ρ (OT file ≈ harmoniser output)
        if rho > 0.95:
            verdict = "✅ CONSISTENT — OT pre-harmonised file is equivalent to harmonizer.sh"
        elif rho > 0.80:
            verdict = "⚠️  MOSTLY CONSISTENT — minor differences"
        else:
            verdict = "❌ DISCORDANT — harmonisation changes results significantly"
    else:
        # Non-harmonised studies: does running the harmoniser actually change things?
        if rho > 0.95:
            verdict = "✅ CONSISTENT — harmonisation has little effect (raw file is fine)"
        elif rho > 0.80:
            verdict = "⚠️  MOSTLY CONSISTENT — minor effect of harmonisation"
        else:
            verdict = "⚠️  DISCORDANT — harmonisation matters for un-harmonised studies"

    print(f"\n  Verdict: {verdict}")
    return {"study_id": study_id, "group": group, "spearman_rho": rho, "p": p,
            "sig_a": (r_a["FDR"] < 0.05).sum(), "sig_b": (r_b["FDR"] < 0.05).sum(),
            "verdict": verdict}


python27 = get_python27()
print(f"\nPython 2.7 for LDSC: {python27}\n")

summary = []
skipped = []

for _, row in sample.iterrows():
    study_id = row["studyId"]
    project  = row["projectId"]
    group    = row["group"]
    n        = int(row["nSamples"]) if pd.notna(row["nSamples"]) else 100_000

    print(f"\n{'#'*60}")
    print(f"  {study_id}  |  Group {group}  |  {row['trait']}")
    print(f"  project={project}  N={n:,}  cases={row['nCases']}")
    print(f"{'#'*60}")

    study_dir = os.path.join(OUT_DIR, study_id)
    os.makedirs(study_dir, exist_ok=True)

    # ── Obtain files ──────────────────────────────────────────────────────────
    hm_path  = None
    raw_path = None

    if group == "A":
        if project == "FINNGEN_R12":
            # Use summarystatsLocation from the parquet (GCS path) instead of EBI FTP
            gcs_path = row.get("summarystatsLocation")
            if not gcs_path:
                print("  No summarystatsLocation for this FinnGen study — skipping.")
                skipped.append({"study_id": study_id, "group": group, "reason": "No GCS path"})
                continue
            fname    = os.path.basename(gcs_path)
            hm_path  = os.path.join(study_dir, fname)
            raw_path = hm_path   # same file used as Path B input
            if not download_finngen(gcs_path, hm_path):
                skipped.append({"study_id": study_id, "group": group, "reason": "FinnGen gsutil failed"})
                continue
        else:
            # GCST Group A: harmonised URL from HM_LOOKUP, raw from one FTP request
            hm_url = HM_LOOKUP.get(study_id)
            if not hm_url:
                print(f"  {study_id} not in harmonised_list.txt — skipping.")
                skipped.append({"study_id": study_id, "group": group, "reason": "Not in HM list"})
                continue
            print(f"  HM URL : {hm_url}")
            raw_url = find_raw_url(study_id)
            print(f"  Raw URL: {raw_url}")
            if not raw_url:
                skipped.append({"study_id": study_id, "group": group, "reason": "No raw URL on FTP"})
                continue
            hm_path  = os.path.join(study_dir, f"{study_id}.h.tsv.gz")
            raw_path = os.path.join(study_dir, f"{study_id}_raw_original.gz")
            if not download_file(hm_url, hm_path):
                skipped.append({"study_id": study_id, "group": group, "reason": "HM download failed"})
                continue
            if not download_file(raw_url, raw_path):
                skipped.append({"study_id": study_id, "group": group, "reason": "Raw download failed"})
                continue

    else:  # Group B — no harmonised file; only raw
        raw_url = find_raw_url(study_id)
        print(f"  Raw URL: {raw_url}")
        if not raw_url:
            skipped.append({"study_id": study_id, "group": group, "reason": "No raw URL on FTP"})
            continue
        raw_path = os.path.join(study_dir, f"{study_id}_raw_original.gz")
        if not download_file(raw_url, raw_path):
            skipped.append({"study_id": study_id, "group": group, "reason": "Raw download failed"})
            continue

    # ── Reformat raw file (used by Path B for all groups; also Path A for Group B) ──
    print("\n  Reformatting raw file ...")
    reformatted_path = reformat_gwas(raw_path, study_dir)
    if reformatted_path is None:
        skipped.append({"study_id": study_id, "group": group, "reason": "Reformat failed"})
        continue

    # ── Path A ────────────────────────────────────────────────────────────────
    if group == "A":
        print("\n  [Path A] Building sumstats from OT/FinnGen harmonised file ...")
        ss_a = build_sumstats_ot(hm_path, n)
    else:
        # Group B: raw file used directly with no harmonisation
        print("\n  [Path A] Building sumstats from raw file (no harmonisation) ...")
        ss_a = build_sumstats_raw(reformatted_path, n)

    if ss_a is None or len(ss_a) < 1000:
        print(f"  Too few SNPs Path A ({len(ss_a) if ss_a is not None else 0}) — skipping.")
        skipped.append({"study_id": study_id, "group": group, "reason": "Too few SNPs (Path A)"})
        continue
    ss_a_path = os.path.join(study_dir, f"{study_id}_path_a.sumstats.gz")
    ss_a.to_csv(ss_a_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(ss_a):,} SNPs → {os.path.basename(ss_a_path)}")

    # ── Path B: raw → harmonizer.sh (same for both groups) ───────────────────
    print("\n  [Path B] Running raw file through harmonizer.sh ...")
    harm_out_dir = run_harmonizer(reformatted_path, study_dir)
    if harm_out_dir is None:
        skipped.append({"study_id": study_id, "group": group, "reason": "Harmonizer failed"})
        continue

    ss_b = harmonized_to_sumstats(harm_out_dir, reformatted_path, n)
    if ss_b is None or len(ss_b) < 1000:
        print(f"  Too few SNPs Path B ({len(ss_b) if ss_b is not None else 0}) — skipping.")
        skipped.append({"study_id": study_id, "group": group, "reason": "Too few SNPs (Path B)"})
        continue
    ss_b_path = os.path.join(study_dir, f"{study_id}_path_b.sumstats.gz")
    ss_b.to_csv(ss_b_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(ss_b):,} SNPs → {os.path.basename(ss_b_path)}")

    # Z-score sanity check
    shared = set(ss_a["SNP"]) & set(ss_b["SNP"])
    if len(shared) > 100:
        z_a = ss_a.set_index("SNP").loc[list(shared), "Z"]
        z_b = ss_b.set_index("SNP").loc[list(shared), "Z"]
        r, _ = stats.pearsonr(z_a, z_b)
        print(f"\n  Z-score Pearson r (shared SNPs, n={len(shared):,}): {r:.4f}")
    else:
        print(f"  WARNING: only {len(shared)} shared SNPs between Path A and B")
        r = None

    # ── LDSC h2-cts ───────────────────────────────────────────────────────────
    print("\n  Running LDSC h2-cts (Path A) ...")
    res_a = run_ldsc(ss_a_path, os.path.join(study_dir, f"{study_id}_ldsc_path_a"), python27)

    print("  Running LDSC h2-cts (Path B) ...")
    res_b = run_ldsc(ss_b_path, os.path.join(study_dir, f"{study_id}_ldsc_path_b"), python27)

    if res_a and res_b:
        result = compare_results(study_id, res_a, res_b, group)
        result["z_pearson_r"] = r
        summary.append(result)
    else:
        skipped.append({"study_id": study_id, "group": group, "reason": "LDSC failed"})


# ── Final summary ─────────────────────────────────────────────────────────
print(f"\n\n{'='*60}")
print("  FINAL SUMMARY")
print(f"{'='*60}")

if summary:
    s = pd.DataFrame(summary)
    print(s[["study_id", "group", "z_pearson_r", "spearman_rho",
              "sig_a", "sig_b", "verdict"]].to_string(index=False))

    for grp, label in [
        ("A", "Group A — pre-harmonised    (expect ρ > 0.95)"),
        ("B", "Group B — NOT harmonised    (expect ρ < 0.95)"),
    ]:
        sub = s[s["group"] == grp]
        if len(sub):
            print(f"\n  {label}:")
            print(f"    Studies analysed      : {len(sub)}")
            print(f"    Consistent (ρ > 0.95) : {(sub['spearman_rho'] > 0.95).sum()} / {len(sub)}")
            print(f"    Mean ρ                : {sub['spearman_rho'].mean():.4f}")

if skipped:
    print(f"\n  Skipped ({len(skipped)}):")
    for sk in skipped:
        print(f"    [Group {sk['group']}] {sk['study_id']}: {sk['reason']}")

print("\nDone.")
