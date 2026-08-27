# -*- coding: utf-8 -*-

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import urllib.request
    from scipy.stats import binomtest
    from statsmodels.stats.multitest import fdrcorrection
    import os
    import re
    import requests
    import time
    import subprocess
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import json
    from scipy.stats import binomtest,mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection
    import glob
    import multiprocessing
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import google.generativeai as genai 
    from dotenv import load_dotenv, find_dotenv
    import zipfile
    import ollama
    from tqdm import tqdm
    import ssl 
    import torch
    import transformers
    import gzip    
    import shutil  
    from transformers import AutoTokenizer, AutoModel
    import torch.nn.functional as F
    import mygene
    from pyfaidx import Fasta
    from scipy.spatial.distance import cosine
    from pyliftover import LiftOver
    import pyranges as pr

    return (
        Fasta,
        LiftOver,
        Path,
        ThreadPoolExecutor,
        as_completed,
        binomtest,
        cosine,
        genai,
        glob,
        gzip,
        json,
        mo,
        multiprocessing,
        np,
        ollama,
        os,
        pd,
        pr,
        re,
        requests,
        shutil,
        ssl,
        subprocess,
        time,
        torch,
        tqdm,
        transformers,
        urllib,
        zipfile,
    )


@app.cell
def _(torch):
    if torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"
    return (device,)


@app.cell
def _(mo):
    mo.md("""
    This notebook reproduces the LDSC cell-type-specific heritability analysis
    from *A single-cell atlas of chromatin accessibility in the human genome* (Zhang et al. 2021).

    All analyses are performed using **GRCh38 / hg38** coordinates.
    """)
    return


@app.cell
def _(mo):
    S3_BASE          = "s3://rejuve-bio/hypothesis-generation-demo"
    GWAS_INPUT_FILE  = mo.ui.text(
        value="C:/Users/Edil/Desktop/hypothesis-generation_demo/data/gwas/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz",
        label="GWAS input file path", 
        full_width=True,
    )
    W_HM3_SNPLIST   = f"{S3_BASE}/ldsc/data/w_hm3.snplist"
    HM3_NO_MHC_LIST = f"{S3_BASE}/data/reference/hm3_no_MHC.list.txt"
    CATLAS_DIR       = "C:/Users/Edil/Desktop/hypothesis-generation_demo/data/catlas_beds"
    CATLAS_URL       = "http://catlas.org/humanenhancer/data/cCREs/"

    mo.vstack([
        mo.md("### Configuration"),
        GWAS_INPUT_FILE,
        mo.md(f"w_hm3.snplist: `{W_HM3_SNPLIST}`"),
    ])
    return (
        CATLAS_DIR,
        CATLAS_URL,
        GWAS_INPUT_FILE,
        HM3_NO_MHC_LIST,
        S3_BASE,
        W_HM3_SNPLIST,
    )


@app.cell
def _(GWAS_INPUT_FILE, os, re):
    _path = GWAS_INPUT_FILE.value
    _basename = os.path.basename(_path)
    _no_ext = _basename
    for _ext in [".tsv.gz", ".txt.gz", ".gz", ".tsv", ".txt", ".csv", ".bgz"]:
        if _no_ext.endswith(_ext):
            _no_ext = _no_ext[: -len(_ext)]
            break
    GWAS_STEM      = re.sub(r"[^A-Za-z0-9_\-]", "_", _no_ext)
    GWAS_FILE      = _path
    SUMSTATS_FILE  = f"data/ldsc_input/{GWAS_STEM}.sumstats.gz"   
    CTS_FILE       = f"data/{GWAS_STEM}_cell_types.cts"
    RESULTS_PREFIX = f"new_results/{GWAS_STEM}_CellTypeSpecific_baseline_v1_weights"

    print(f"GWAS stem      : {GWAS_STEM}")
    print(f"GWAS file      : {GWAS_FILE}")
    print(f"Sumstats file  : {SUMSTATS_FILE}")   
    print(f"Results prefix : {RESULTS_PREFIX}")
    return CTS_FILE, GWAS_FILE, GWAS_STEM, RESULTS_PREFIX, SUMSTATS_FILE


@app.cell
def _(GWAS_FILE, os, pd):
    def load_and_standardize(file_path):

        _filename = os.path.basename(file_path)
        rename_map = {

            'hm_rsid': 'rsid', 'hm_variant_id': 'rsid', 'rs_dbSNP147': 'rsid', 'SNP': 'rsid',

            'hm_chrom': 'chr', 'chromosome': 'chr', 'CHR': 'chr',

            'hm_pos': 'pos', 'base_pair_location': 'pos', 'POS_GRCh37': 'pos',

            'p_value': 'p', 'P-value': 'p', 'pval': 'p', 'Pvalue': 'p',

            'hm_beta': 'beta', 'beta': 'beta', 'Effect_A2': 'beta',

            'Freq_A2': 'maf', 'minor_AF': 'maf', 'Freq': 'maf', 'af': 'maf',



            'A1': 'ref', 'A2': 'alt', 

            'hm_other_allele': 'ref', 'hm_effect_allele': 'alt',

            'other_allele': 'ref', 'effect_allele': 'alt'

        }


        _chunks = []



        _reader = pd.read_csv(file_path, sep=r'\s+', engine='c', compression='infer', chunksize=250000)


        for _chunk in _reader:



            _existing = {k: v for k, v in rename_map.items() if k in _chunk.columns}

            _chunk = _chunk.rename(columns=_existing)





            _chunk = _chunk.loc[:, ~_chunk.columns.duplicated(keep='last')]




            _needed = ['rsid', 'chr', 'pos', 'p', 'ref', 'alt', 'beta', 'maf']

            _cols_present = [c for c in _needed if c in _chunk.columns]

            _chunk = _chunk[_cols_present]




            if 'p' in _chunk.columns:

                _chunk['p'] = pd.to_numeric(_chunk['p'], errors='coerce')

                _chunk = _chunk[_chunk['p'] < 5e-8].dropna(subset=['rsid', 'p']).copy()




            if 'maf' in _chunk.columns:

                _chunk['maf'] = pd.to_numeric(_chunk['maf'], errors='coerce')

                _chunk = _chunk[_chunk['maf'] >= 0.01].copy()


            if not _chunk.empty:

                _chunks.append(_chunk)


        if not _chunks:

            print(" ERROR: No SNPs passed the MAF 0.01 and P-value filters.")

            return pd.DataFrame()




        df = pd.concat(_chunks, ignore_index=True)

        df = df.sort_values('p', ascending=True).drop_duplicates(subset='rsid', keep='first')



        df['chr'] = pd.to_numeric(df['chr'], errors='coerce').fillna(0).astype(int).astype(str)

        df['pos'] = pd.to_numeric(df['pos'], errors='coerce').fillna(0).astype(int)




        if 'ref' not in df.columns:

             print(" Warning: 'ref' column not found. Checking raw columns...")

             print(f"Available: {list(df.columns)}")


        print(f" SUCCESS: {len(df)} variants identified.")

        return df


    current_gwas = load_and_standardize(GWAS_FILE)
    return (current_gwas,)


@app.cell
def _(mo):
    mo.md("""
    ## Add Biological Infrastructure
    """)
    return


@app.cell
def _(Path, os, ssl, urllib, zipfile):

    ssl._create_default_https_context = ssl._create_unverified_context



    BIN_DIR = Path("bin")

    BIN_DIR.mkdir(exist_ok=True)



    PLINK_BIN = BIN_DIR / "plink.exe"

    CHAIN_FILE = Path("data/reference/hg19ToHg38.over.chain.gz")




    if not PLINK_BIN.exists():

        _plink_url = "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip"

        _zip_path = "plink_windows.zip"

        urllib.request.urlretrieve(_plink_url, _zip_path)





        with zipfile.ZipFile(_zip_path, 'r') as zip_ref:

            zip_ref.extractall(str(BIN_DIR))





        if os.path.exists(_zip_path):

            os.remove(_zip_path)

        print("  Result: PLINK.exe is ready.")




    if not CHAIN_FILE.exists():

        os.makedirs(CHAIN_FILE.parent, exist_ok=True)

        print("Action: Downloading hg19ToHg38 chain file...")

        urllib.request.urlretrieve(

            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", 

            str(CHAIN_FILE)

        )

        print("  Result: Chain file is ready.")

    return CHAIN_FILE, PLINK_BIN


@app.cell
def download_genome_reference(Path, gzip, os, requests, shutil):

    save_dir = Path("data/reference/GRCh38")

    save_dir.mkdir(parents=True, exist_ok=True) 





    GENOME_FASTA_PATH = save_dir / "hg38_analysis_set.fa"

    _gz_file = save_dir / "hg38.analysisSet.fa.gz"


    if not GENOME_FASTA_PATH.exists():



        _genome_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"



        try:



            _response = requests.get(_genome_url, stream=True, timeout=300)

            _response.raise_for_status() 



            with open(_gz_file, 'wb') as f:

                for _chunk in _response.iter_content(chunk_size=1024*1024): # 1MB chunks

                    if _chunk:

                        f.write(_chunk)





            with gzip.open(_gz_file, 'rb') as f_in:

                with open(GENOME_FASTA_PATH, 'wb') as f_out:

                    shutil.copyfileobj(f_in, f_out)





            if os.path.exists(_gz_file):

                os.remove(_gz_file)




        except Exception as e:

            print(f"  CRITICAL ERROR: {e}")

            if os.path.exists(_gz_file):

                os.remove(_gz_file) 

    else:

        print(f"  Result: Reference file '{GENOME_FASTA_PATH.name}' already exists.")
    return (GENOME_FASTA_PATH,)


@app.cell
def _(CATLAS_DIR, Path, mo, os, requests, tqdm):
    _catlas_path = Path(CATLAS_DIR)
    _catlas_path.mkdir(parents=True, exist_ok=True)
    _metadata_filename = "catlas_cell_types.txt"
    _search_paths = [
        _metadata_filename,
        os.path.join("LLM_Integration", _metadata_filename),
        os.path.join("..", "LLM_Integration", _metadata_filename)
    ]

    _valid_meta_path = next((p for p in _search_paths if os.path.exists(p)), None)
    sync_report_ui = None 
    if _valid_meta_path is None:
        sync_report_ui = mo.md(f"###  Error\nMetadata file `{_metadata_filename}` not found on disk.")
    else:
        with open(_valid_meta_path, "r") as _f_in:
            _cell_types = [line.strip() for line in _f_in if line.strip()]
        for _ct_name in tqdm(_cell_types, desc="Syncing CATlas"):
            _local_file = _catlas_path / f"{_ct_name}.bed"
            if not _local_file.exists():
                _url = f"http://catlas.org/humanenhancer/data/cCREs/{_ct_name}.bed"
                try:
                    _res = requests.get(_url, timeout=30)
                    if _res.status_code == 200:
                        with open(_local_file, 'wb') as _data_out:
                            _data_out.write(_res.content)
                except:
                    continue
        _final_count = len([f for f in os.listdir(_catlas_path) if f.endswith(".bed")])
        sync_report_ui = mo.md(
            f"Successfully verified **{_final_count}** cell types in `{_catlas_path.name}`.\n"
        )
    sync_report_ui
    return


@app.cell
def _(Path):

    _eur_dir = Path("C:/Users/Edil/Desktop/hypothesis-generation_demo/data/EUR")

    _prefix = "1000G.EUR.QC"
    _found_chroms = []

    for _ch in range(1, 23):

        if (_eur_dir / f"{_prefix}.{_ch}.bim").exists():

            _found_chroms.append(_ch)





    if len(_found_chroms) == 22:

        EUR_REF_TEMPLATE = str(_eur_dir / _prefix)

    else:

        print(f" WARNING: Found only {len(_found_chroms)} / 22 chromosomes.")

        EUR_REF_TEMPLATE = "INCOMPLETE"
    return (EUR_REF_TEMPLATE,)



@app.cell
def _(mo):
    mo.md("""
    ## Liftover & Allele Reorientation
    """)
    return


@app.cell
def _(CHAIN_FILE, LiftOver, pd):
    def run_liftover(df_hg19):

        lo = LiftOver(str(CHAIN_FILE))

        lifted_data = []


        for _, row in df_hg19.iterrows():



            _c = str(row['chr'])

            _chrom = f"chr{_c}" if not _c.startswith('chr') else _c





            res = lo.convert_coordinate(_chrom, int(row['pos']) - 1)

            if res:



                lifted_data.append({

                    'rsid': row['rsid'],

                    'pos_hg38': res[0][1] + 1 

                })


        if not lifted_data:

            return pd.DataFrame()


        lifted_df = pd.DataFrame(lifted_data)



        df_hg38 = df_hg19.merge(lifted_df, on='rsid', how='inner')

        return df_hg38


    def reorient_alleles(df, genome_reader):

        def _fix_row(row):

            _c = str(row['chr'])

            _chrom = f"chr{_c}" if not _c.startswith('chr') else _c

            _pos = int(row['pos_hg38']) 

            try:



                _genome_seq = genome_reader[_chrom][_pos-1:_pos].seq.upper()

                if _genome_seq == row['ref']:

                    return row['ref'], row['alt'], row['beta']

                elif _genome_seq == row['alt']:



                    return row['alt'], row['ref'], -row['beta']

                else:

                    return None, None, None 

            except:

                return None, None, None


        results = df.apply(_fix_row, axis=1, result_type='expand')

        results.columns = ['ref_fixed', 'alt_fixed', 'beta_fixed']



        df_out = df.copy()

        df_out['ref'] = results['ref_fixed']

        df_out['alt'] = results['alt_fixed']

        df_out['beta'] = results['beta_fixed']




        before = len(df_out)

        df_out = df_out.dropna(subset=['ref', 'alt'])

        print(f"  Result: {len(df_out)} oriented variants ready.")

        return df_out

    return (run_liftover,)


@app.cell
def _(mo):
    mo.md("""
    ## NT-2.5B Scoring Engine
    """)
    return


@app.cell
def _(device, np, torch, transformers):

    _MODEL_ID = "InstaDeepAI/nucleotide-transformer-2.5b-multi-species"

    tokenizer = transformers.AutoTokenizer.from_pretrained(_MODEL_ID, trust_remote_code=True)

    model = transformers.AutoModel.from_pretrained(_MODEL_ID, trust_remote_code=True).to(device).eval()


    def score_variants(seq_df, batch_size=2):

        _scores = []

        _embeddings = []

        print(f"  Action: Scoring {len(seq_df)} variants via NT-2.5B...")


        with torch.no_grad():

            for i in range(0, len(seq_df), batch_size):

                _batch = seq_df.iloc[i : i + batch_size]

                _ref_in = tokenizer(_batch['seq_ref'].tolist(), return_tensors="pt", padding=True).to(device)

                _alt_in = tokenizer(_batch['seq_alt'].tolist(), return_tensors="pt", padding=True).to(device)


                _ref_emb = model(**_ref_in).last_hidden_state.mean(dim=1)

                _alt_emb = model(**_alt_in).last_hidden_state.mean(dim=1)

                _batch_scores = 1 - torch.nn.functional.cosine_similarity(_ref_emb, _alt_emb)

                _scores.extend(_batch_scores.cpu().numpy().tolist())

                _embeddings.append(_ref_emb.cpu().numpy())


        return _scores, np.vstack(_embeddings)

    return (score_variants,)


@app.cell
def _(pd):
    def extract_dna_windows(df_hg38, genome_reader, window_size=1000):

        _half = window_size // 2

        _rows = []

        for _, row in df_hg38.iterrows():

            _c = str(row['chr'])

            _chrom = f"chr{_c}" if not _c.startswith('chr') else _c

            _pos = int(row['pos_hg38'])



            try:

                _full_seq = genome_reader[_chrom][_pos - 1 - _half : _pos + _half].seq.upper()

                if len(_full_seq) < window_size:

                    continue

                _seq_ref = _full_seq[:_half] + str(row['ref']) + _full_seq[_half+1:]

                _seq_alt = _full_seq[:_half] + str(row['alt']) + _full_seq[_half+1:]

                _rows.append({

                    'rsid': row['rsid'],

                    'chr': row['chr'],

                    'pos_hg38': _pos,

                    'seq_ref': _seq_ref,

                    'seq_alt': _seq_alt

                })

            except Exception as e:



                continue 


        _result_df = pd.DataFrame(_rows)

        return _result_df




    return (extract_dna_windows,)


@app.cell
def universal_discovery_engine(
    CATLAS_DIR,
    EUR_REF_TEMPLATE,
    Fasta,
    GENOME_FASTA_PATH,
    GWAS_STEM,
    PLINK_BIN,
    binomtest,
    cosine,
    current_gwas,
    extract_dna_windows,
    mo,
    np,
    os,
    pd,
    pr,
    run_liftover,
    score_variants,
    subprocess,
    tqdm,
):


    _res_dir = "data/results"
    os.makedirs(_res_dir, exist_ok=True)
    _v16_scores_path = f"{_res_dir}/{GWAS_STEM}_FINAL_V16_AI_SCORES.csv"
    _v16_table_path = f"{_res_dir}/{GWAS_STEM}_FINAL_V16_LEADERBOARD.csv"
    _genome_bp = 3_100_000_000

    # --- . PRUNING (hg19) ---
    _all_pruned_hits = []
    for _ch in range(1, 23):
        _bim_p = f"{EUR_REF_TEMPLATE}.{_ch}.bim"
        if not os.path.exists(_bim_p): continue
        _bim_df = pd.read_csv(_bim_p, sep=r"\s+", header=None, usecols=[1], names=['rsid'], engine='python')
        _bim_ids = set(_bim_df['rsid'].astype(str).str.strip())
        _vetted_ch = current_gwas[current_gwas['rsid'].isin(_bim_ids)].copy()
        if not _vetted_ch.empty:
            _ch_rs = f"data/gwas/rs_v16_ch{_ch}.txt"; _vetted_ch['rsid'].to_csv(_ch_rs, index=False, header=False)
            _out = f"data/gwas/p_v16_ch{_ch}"
            subprocess.run([str(PLINK_BIN), "--bfile", f"{EUR_REF_TEMPLATE}.{_ch}", "--extract", _ch_rs, "--indep-pairwise", "1000", "100", "0.1", "--out", _out], capture_output=True)
            if os.path.exists(f"{_out}.prune.in"):
                with open(f"{_out}.prune.in", 'r') as _f_win: _win_rsids = _f_win.read().splitlines()
                _all_pruned_hits.append(_vetted_ch[_vetted_ch['rsid'].isin(_win_rsids)])

    _pruned_hg19 = pd.concat(_all_pruned_hits, ignore_index=True)
    _snp_pr_hg19 = pr.PyRanges(pd.DataFrame({'Chromosome': _pruned_hg19['chr'].apply(lambda x: f"chr{x}" if not str(x).startswith('chr') else str(x)), 'Start': _pruned_hg19['pos'] - 1, 'End': _pruned_hg19['pos'], 'rsid': _pruned_hg19['rsid']}))
    _all_hits_accum = []
    _beds = [f for f in os.listdir(CATLAS_DIR) if f.endswith(".bed")]
    for _f_name in tqdm(_beds, desc="Checking Specificity"):
        _ct_inner = _f_name.replace(".bed", "")
        _ct_bed = pd.read_csv(os.path.join(CATLAS_DIR, _f_name), sep='\t', header=None, usecols=[0,1,2], names=['Chromosome', 'Start', 'End'])
        _overlaps = _snp_pr_hg19.overlap(pr.PyRanges(_ct_bed)).as_df()
        if not _overlaps.empty:
            for _r in _overlaps['rsid'].unique(): _all_hits_accum.append({'rsid': _r, 'Cell_Type': _ct_inner})

    _breadth = pd.DataFrame(_all_hits_accum)['rsid'].value_counts()
    _promisc_ids = _breadth[_breadth > (len(_beds) * 0.20)].index 
    _df_v16_filtered = _pruned_hg19[~_pruned_hg19['rsid'].isin(_promisc_ids)].copy()

    # --- . TRANSLATE & REORIENT ---
    _df_h38 = run_liftover(_df_v16_filtered)
    _reader = Fasta(str(GENOME_FASTA_PATH))
    _comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    def _fix_strand(r):
        try:
            _g = _reader[f"chr{r['chr']}"][int(r['pos_hg38'])-1:int(r['pos_hg38'])].seq.upper()
            if _g == r['ref']: return r['ref'], r['alt'], r['beta']
            if _g == r['alt']: return r['alt'], r['ref'], -r['beta']
            if _g == _comp.get(r['ref']): return _g, _comp.get(r['alt']), r['beta']
            return r['ref'], r['alt'], r['beta']
        except: return r['ref'], r['alt'], r['beta']

    _res_al = _df_h38.apply(_fix_strand, axis=1, result_type='expand')
    _res_al.columns = ['ref_fix', 'alt_fix', 'beta_fix']
    _df_ready = _df_h38.assign(ref=_res_al['ref_fix'], alt=_res_al['alt_fix'], beta=_res_al['beta_fix'])


    _scored_list = []
    if os.path.exists(_v16_scores_path) and os.path.getsize(_v16_scores_path) > 100:
        _scored_list = pd.read_csv(_v16_scores_path).to_dict('records')

    if len(_scored_list) < len(_df_ready):
        _start = len(_scored_list)
        _seq_df = extract_dna_windows(_df_ready.iloc[_start:], _reader)
        for _i in range(len(_seq_df)):
            _row_in = _seq_df.iloc[_i:_i+1]
            _d_val, _e_val = score_variants(_row_in, batch_size=1)
            _scored_list.append({'rsid': _row_in['rsid'].values[0], 'impact_score': _d_val[0], 'embedding_json': str(list(_e_val[0]))})
            if _i % 10 == 0: pd.DataFrame(_scored_list).to_csv(_v16_scores_path, index=False)

    _impact_v16 = pd.DataFrame(_scored_list).drop_duplicates('rsid')
    _results = _df_ready.merge(_impact_v16, on='rsid', how='inner')
    _results['vec_obj'] = _results['embedding_json'].apply(lambda x: np.array(eval(x)) if isinstance(x, str) else np.array(x))
    _N_total = len(_results)
    _snp_pr = pr.PyRanges(pd.DataFrame({'Chromosome': _results['chr'].apply(lambda x: f"chr{x}" if not str(x).startswith('chr') else str(x)), 'Start': _results['pos_hg38'] - 1, 'End': _results['pos_hg38'], 'rsid': _results['rsid']}))
    _leaderboard = []
    for _f_name in tqdm(_beds, desc="Final Scanning"):
        _ct_bed_df = pd.read_csv(os.path.join(CATLAS_DIR, _f_name), sep='\t', header=None, usecols=[0,1,2], names=['Chromosome', 'Start', 'End'])
        _hits_overlaps = _snp_pr.overlap(pr.PyRanges(_ct_bed_df)).as_df()
        if not _hits_overlaps.empty:
            _prob_bg = int((pr.PyRanges(_ct_bed_df).End - pr.PyRanges(_ct_bed_df).Start).sum()) / _genome_bp
            _test = binomtest(len(_hits_overlaps), _N_total, _prob_bg, alternative='greater')
            # Peak In vs Out Contrast
            _in_m = _results['rsid'].isin(set(_hits_overlaps['rsid']))
            _mean_in = np.mean(_results.loc[_in_m, 'vec_obj'].tolist(), axis=0)
            _mean_out = np.mean(_results.loc[~_in_m, 'vec_obj'].tolist(), axis=0)
            _dist = cosine(_mean_in, _mean_out) if _in_m.any() and (~_in_m).any() else 0
            _leaderboard.append({'Cell_Type': _f_name.replace(".bed",""), 'Hits': len(_hits_overlaps), 'P_Value': round(_test.pvalue, 6), 'AI_Validation': round(_dist, 4)})

    _df_final = pd.DataFrame(_leaderboard).sort_values('P_Value').reset_index(drop=True)
    _df_final.to_csv(_v16_table_path, index=False)
    mo.vstack([
        mo.md(f"Discovery Leaderboard: {GWAS_STEM}"),
        mo.ui.table(_df_final.head(15))
    ])
    return


@app.cell
def _(mo):
    mo.md("""
    ## 0. Setup: download and configure LDSC
    """)
    return


@app.cell
def _(Path, json, os, subprocess):
    TOOLS_DIR = Path("tools")
    LDSC_DIR  = TOOLS_DIR / "ldsc"
    TOOLS_DIR.mkdir(exist_ok=True)

    _env_check = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if "ldsc27" not in _env_check.stdout:
        subprocess.run(["conda", "create", "-n", "ldsc27", "python=2.7", "-y"], check=True)

    _conda_json = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    _envs = json.loads(_conda_json.stdout)["envs"]
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
    python27_path = os.path.join(ldsc27_path, "bin", "python")

    print(f"LDSC environment ready")
    print(f"Python 2.7 : {python27_path}")
    print(f"LDSC dir   : {LDSC_DIR}")
    return ldsc27_path, python27_path


@app.cell
def _(mo):
    mo.md("""
    ## 1. Discover cell types and resolve BED sources
    """)
    return


@app.cell
def _(CATLAS_DIR, CATLAS_URL, S3_BASE, os, pd, re, subprocess, urllib):
    for _d in ["data/peaks", "data/reference", "data/gwas", "data/beds"]:
        os.makedirs(_d, exist_ok=True)

    BED_SEARCH_DIRS = [CATLAS_DIR, "data/beds"]

    def _sanitize(name):
        return re.sub(r"[^A-Za-z0-9_\-]", "_", name)

    def _ensure_catlas():
        if os.path.exists(CATLAS_DIR):
            _beds = [f for f in os.listdir(CATLAS_DIR) if f.endswith(".bed")]
            if _beds:
                print(f"{CATLAS_DIR} has {len(_beds)} BED files, skipping download")
                return
        os.makedirs(CATLAS_DIR, exist_ok=True)
        subprocess.run(
            ["wget", "-r", "-np", "-nH", "--cut-dirs=3",
             "-R", "index.html*", "-P", CATLAS_DIR, CATLAS_URL],
            check=True,
        )
        print(f"Downloaded {len([f for f in os.listdir(CATLAS_DIR) if f.endswith('.bed')])} BED files")

    def _find_local_bed(name):
        for d in BED_SEARCH_DIRS:
            p = os.path.join(d, f"{name}.bed")
            if os.path.exists(p):
                return p
        return None

    def _peak_to_bed(ct):
        peak_txt = f"{S3_BASE}/data/peaks/{ct}.peak.annotation.txt"
        bed_out  = f"data/beds/{ct}.bed"
        if os.path.exists(bed_out):
            return bed_out
        _peaks = pd.read_csv(peak_txt, sep="\t")
        _peaks[["seqnames", "start", "end"]].rename(
            columns={"seqnames": "chr"}
        ).to_csv(bed_out, sep="\t", index=False, header=False)
        return bed_out

    _ensure_catlas()

    raw_names = set()
    for _sd in BED_SEARCH_DIRS:
        if os.path.exists(_sd):
            for _f in os.listdir(_sd):
                if _f.endswith(".bed") and not _f.startswith("."):
                    raw_names.add(os.path.splitext(_f)[0])

    cell_type_beds = {}
    _seen = {}
    for _raw in sorted(raw_names):
        _safe = _sanitize(_raw)
        if _safe in _seen:
            _safe += "_2"
        _seen[_safe] = _raw
        _local = _find_local_bed(_raw)
        cell_type_beds[_safe] = _local if _local else _peak_to_bed(_raw)

    all_cell_types = sorted(cell_type_beds.keys())
    print(f"{len(all_cell_types)} cell types ready")

    for _url, _dst in [
        (
            "https://zenodo.org/records/10515792/files/GRCh38.tgz?download=1",
            "data/reference/GRCh38.tgz",
        ),
        (
            "https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1",
            "data/reference/hm3_no_MHC.list.txt",
        ),
    ]:
        if not os.path.exists(_dst):
            print(f"Downloading {_dst}...")
            urllib.request.urlretrieve(_url, _dst)

    print("All sources resolved")
    return all_cell_types, cell_type_beds


@app.cell
def _(mo):
    mo.md("""
    ## 2. Extract reference LD panels
    """)
    return


@app.cell
def _(S3_BASE, os, subprocess):
    if not os.path.exists("data/reference/GRCh38"):
        subprocess.run(
            ["tar", "-xzf", "data/reference/GRCh38.tgz", "-C", "data/reference"],
            check=True,
        )

    for _tgz, _base, _check in [
        (
            f"{S3_BASE}/data/reference/GRCh38/plink_files.tgz",
            "data/reference/GRCh38",
            "plink_files/1000G.EUR.hg38.1.bim",
        ),
        (
            f"{S3_BASE}/data/reference/GRCh38/weights.tgz",
            "data/reference/GRCh38",
            "weights",
        ),
    ]:
        if not os.path.exists(os.path.join(_base, _check)) and os.path.exists(_tgz):
            subprocess.run(["tar", "-xzf", _tgz, "-C", _base], check=True)

    _critical = "data/reference/GRCh38/plink_files/1000G.EUR.hg38.1.bim"
    print("Reference files ready" if os.path.exists(_critical) else f"ERROR: missing {_critical}")
    return


@app.cell
def _(mo):
    mo.md("""
    ## 3. Reformat GWAS to harmonizer-compatible format
    """)
    return


@app.cell
def _(GWAS_FILE, SUMSTATS_FILE, glob, os, pd):
    import gzip as _gzip
    import shutil as _shutil

    _lower   = GWAS_FILE.lower()
    _is_gz   = _lower.endswith(".gz") or _lower.endswith(".bgz")
    _comp    = "gzip" if _is_gz else None
    _open_fn = _gzip.open if _is_gz else open

    _skip = 0
    _peek_line = ""
    with _open_fn(GWAS_FILE, "rt") as _fh:
        for _line in _fh:
            if not _line.startswith("##"):
                _peek_line = _line
                break
    _sep = "\t" if "\t" in _peek_line else " "

    os.makedirs("data/gwas_reformatted", exist_ok=True)
    REFORMATTED_GWAS = f"data/gwas_reformatted/{os.path.basename(GWAS_FILE)}.gz"

    if os.path.exists(SUMSTATS_FILE):
        print(f"Sumstats already exists, skipping reformat: {SUMSTATS_FILE}")
    elif os.path.exists(REFORMATTED_GWAS):
        print(f"Reformatted file already exists: {REFORMATTED_GWAS}")
    else:
        _df = pd.read_csv(GWAS_FILE, sep=_sep, compression=_comp, skiprows=_skip, low_memory=False)
        _df.columns = [c.lstrip("#").strip() for c in _df.columns]
        print(f"  Columns: {list(_df.columns)}")

        _cl_check = {c.lower() for c in _df.columns}
        if "z" in _cl_check and "snp" in _cl_check and "beta" not in _cl_check and "logor" not in _cl_check:
            print("  File is already in LDSC format — copying directly to sumstats")
            os.makedirs("data/ldsc_input", exist_ok=True)
            _shutil.copy(GWAS_FILE, SUMSTATS_FILE)
            print(f"  Copied to {SUMSTATS_FILE}")
            REFORMATTED_GWAS = GWAS_FILE
        else:
            _cl = {c.lower(): c for c in _df.columns}

            _snp_col  = _cl.get("markername") or _cl.get("snptestid") or _cl.get("id") or _cl.get("snp") or _cl.get("rsid") or _cl.get("variant_id")
            _chr_col  = _cl.get("chr") or _cl.get("chrom") or _cl.get("chromosome") or _cl.get("#chrom") or _cl.get("hm_chrom")
            _pos_col  = _cl.get("pos") or _cl.get("bp") or _cl.get("position") or _cl.get("bp_hg19") or _cl.get("base_pair_location")
            _a1_col   = _cl.get("a1") or _cl.get("effect_allele") or _cl.get("alt") or _cl.get("hm_effect_allele")
            _a2_col   = _cl.get("a2") or _cl.get("noneffect_allele") or _cl.get("other_allele") or _cl.get("ref") or _cl.get("hm_other_allele")
            _beta_col = _cl.get("logor") or _cl.get("log_or") or _cl.get("beta") or _cl.get("b") or _cl.get("hm_beta")
            _se_col   = _cl.get("se_gc") or _cl.get("stderrlogor") or _cl.get("se") or _cl.get("stderr") or _cl.get("standard_error")
            _p_col    = _cl.get("p-value_gc") or _cl.get("pvalue") or _cl.get("p_value") or _cl.get("p-value") or _cl.get("p")
            _n_col    = _cl.get("n_samples") or _cl.get("neff") or _cl.get("n") or _cl.get("n_total")

            _rename = {}
            if _snp_col:  _rename[_snp_col]  = "snp"
            if _a1_col:   _rename[_a1_col]   = "a1"
            if _a2_col:   _rename[_a2_col]   = "a2"
            if _beta_col: _rename[_beta_col] = "beta"
            if _se_col:   _rename[_se_col]   = "se"
            if _p_col:    _rename[_p_col]    = "p"
            if _n_col:    _rename[_n_col]    = "n"
            if _chr_col:  _rename[_chr_col]  = "chr"
            if _pos_col:  _rename[_pos_col]  = "pos"
            _df = _df.rename(columns=_rename)

            if "chr" not in _df.columns or "pos" not in _df.columns:
                print("  No chr/pos columns — looking up from bim files via rs IDs...")
                _bim_map = pd.concat([
                    pd.read_csv(f, sep="\t", header=None, names=["chr","snp","cm","pos","a1b","a2b"])[["snp","chr","pos"]]
                    for f in sorted(glob.glob("data/reference/GRCh38/plink_files/1000G.EUR.hg38.*.bim"))
                ]).drop_duplicates("snp").set_index("snp")
                _df["chr"] = _df["snp"].map(_bim_map["chr"])
                _df["pos"] = _df["snp"].map(_bim_map["pos"])
                _df = _df.dropna(subset=["chr","pos"])
                _df["chr"] = _df["chr"].astype(int).astype(str)
                _df["pos"] = _df["pos"].astype(int).astype(str)
                print(f"  Mapped {len(_df):,} variants with chr/pos")

            if "chr" in _df.columns:
                _df["chr"] = _df["chr"].astype(str).str.replace("chr","",regex=False)

            _keep = [c for c in ["chr","pos","snp","a1","a2","beta","se","p","n"] if c in _df.columns]
            _dropna_cols = [c for c in ["a1","a2","beta","p"] if c in _df.columns]
            _df[_keep].dropna(subset=_dropna_cols).to_csv(
                REFORMATTED_GWAS, sep="\t", index=False, compression="gzip"
            )
            print(f"  Written {len(_df):,} variants to {REFORMATTED_GWAS}")
    return (REFORMATTED_GWAS,)


@app.cell
def _(mo):
    mo.md("""
    ## 4. Setup harmonization workflow
    """)
    return


@app.cell
def _(Path, os):
    HARMONIZER_CODE_REPO = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser"
    HARMONIZER_REF_DIR   = "/mnt/hdd_1/abdu/gwas-sumstats-harmoniser/data/gwas_harm_ref"
    harmonizer_script    = Path(HARMONIZER_CODE_REPO) / "harmonizer.sh"

    nextflow_env = {
        **os.environ,
        "PATH": "/mnt/hdd_1/rediet/hypothesis-generation-demo:/mnt/hdd_1/rediet/jdk-17/bin:" + os.environ.get("PATH", ""),
        "JAVA_HOME": "/mnt/hdd_1/rediet/jdk-17",
    }

    if not harmonizer_script.exists():
        print(f"WARNING: harmonizer.sh not found at {harmonizer_script}")
        harmonizer_ready = False
    elif not os.path.isdir(HARMONIZER_REF_DIR):
        print(f"WARNING: Reference directory not found at {HARMONIZER_REF_DIR}")
        harmonizer_ready = False
    else:
        print(f"Harmonizer configuration found")
        print(f"  Script   : {harmonizer_script}")
        print(f"  Reference: {HARMONIZER_REF_DIR}")
        harmonizer_ready = True
    return (
        HARMONIZER_CODE_REPO,
        HARMONIZER_REF_DIR,
        harmonizer_ready,
        harmonizer_script,
        nextflow_env,
    )


@app.cell
def _(mo):
    mo.md("""
    ## 5. Harmonize GWAS summary statistics
    """)
    return


@app.cell
def _(
    HARMONIZER_CODE_REPO,
    HARMONIZER_REF_DIR,
    REFORMATTED_GWAS,
    harmonizer_ready,
    harmonizer_script,
    nextflow_env,
    os,
    subprocess,
):
    os.makedirs("data/harmonized", exist_ok=True)

    harmonized_found = False
    harmonized_output_dir = None

    if os.path.exists("data/harmonized"):
        for item in os.listdir("data/harmonized"):
            if os.path.isdir(os.path.join("data/harmonized", item)):
                _final_dir = os.path.join("data/harmonized", item, "final")
                if os.path.exists(_final_dir):
                    harmonized_found = True
                    harmonized_output_dir = os.path.join("data/harmonized", item)
                    break

    _ssf = os.path.join(
        os.path.dirname(REFORMATTED_GWAS),
        os.path.basename(REFORMATTED_GWAS).replace(".gz", ".tsv.gz")
    )
    if not harmonized_found and os.path.exists(_ssf):
        print(f"SSF file exists, skipping harmonization: {_ssf}")
        harmonized_found = True
        harmonized_output_dir = "data/harmonized"

    if harmonized_found:
        print(f"Harmonized GWAS file already exists: {harmonized_output_dir}")
    elif not harmonizer_ready:
        print("Skipping harmonization - harmonizer not configured")
        harmonized_output_dir = None
    else:
        os.makedirs("data/harmonized", exist_ok=True)
        log_file = os.path.abspath(f"data/harmonized/harmonizer_{os.path.basename(REFORMATTED_GWAS)}.log")
        input_abs = os.path.abspath(REFORMATTED_GWAS)

        os.chdir("data/harmonized")

        try:
            with open(log_file, "w") as _log:
                _result = subprocess.run([
                    "bash",
                    str(harmonizer_script),
                    "--input",     input_abs,
                    "--build",     "GRCh38",
                    "--ref",       HARMONIZER_REF_DIR,
                    "--code-repo", HARMONIZER_CODE_REPO,
                    "--threshold", "0.99"
                ], check=False, text=True, env=nextflow_env,
                   stdout=_log, stderr=_log)

            print(f"Return code: {_result.returncode}")
            print(f"Log: {log_file}")

            harmonized_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
            if harmonized_dirs:
                harmonized_output_dir = os.path.abspath(max(harmonized_dirs))
                print(f"Harmonization complete: {harmonized_output_dir}")
            else:
                harmonized_output_dir = None

        except subprocess.CalledProcessError as e:
            print(f"ERROR: Harmonization failed")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            harmonized_output_dir = None
        finally:
            os.chdir("../..")
    return (harmonized_output_dir,)


@app.cell
def _(mo):
    mo.md("""
    ## 6. Convert harmonized output to LDSC format
    """)
    return


@app.cell
def _(
    REFORMATTED_GWAS,
    SUMSTATS_FILE,
    W_HM3_SNPLIST,
    harmonized_output_dir,
    os,
    pd,
):
    import yaml as _yaml

    os.makedirs("data/ldsc_input", exist_ok=True)

    def _find_n_deep(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                low_key = str(k).lower()
                if (('sample' in low_key and 'size' in low_key) or low_key == 'n') and isinstance(v, (int, float)):
                    return int(v)
            for v in obj.values():
                res = _find_n_deep(v)
                if res: return res
        elif isinstance(obj, list):
            for item in obj:
                res = _find_n_deep(item)
                if res: return res
        return None

    if harmonized_output_dir is None:
        print("Skipping - no harmonized data available")
    elif os.path.exists(SUMSTATS_FILE):
        print(f"LDSC sumstats already exists, skipping: {SUMSTATS_FILE}")
    else:
        _final_dir = os.path.join(harmonized_output_dir, "final")
        _files = [f for f in os.listdir(_final_dir) if f.endswith(".tsv.gz")] if os.path.exists(_final_dir) else []

        if not _files:
            print("No final dir — using SSF file directly...")
            _ssf = os.path.join(
                os.path.dirname(REFORMATTED_GWAS),
                os.path.basename(REFORMATTED_GWAS).replace(".gz", ".tsv.gz")
            )
            if os.path.exists(_ssf):
                _files = [os.path.basename(_ssf)]
                _final_dir = os.path.dirname(_ssf)
            else:
                print(f"WARNING: SSF file not found: {_ssf}")

        if not _files:
            print("WARNING: No harmonized or SSF file found")
        else:
            _harmonized_file = os.path.join(_final_dir, _files[0])
            print(f"Converting {_harmonized_file} to LDSC format...")

            _df = pd.read_csv(_harmonized_file, sep="\t", compression="gzip")
            print(f"  Columns: {list(_df.columns)}")

            _df = _df.rename(columns={
                "rsid":            "SNP",
                "effect_allele":   "A1",
                "other_allele":    "A2",
                "beta":            "BETA",
                "standard_error":  "SE",
                "p_value":         "P",
            })

            _df["BETA"] = pd.to_numeric(_df["BETA"], errors="coerce")
            _df["SE"]   = pd.to_numeric(_df["SE"],   errors="coerce")
            _df["A1"]   = _df["A1"].str.upper()
            _df["A2"]   = _df["A2"].str.upper()
            _df["Z"]    = _df["BETA"] / _df["SE"]

            if "N" not in _df.columns:
                _detected_n = None
                for _col in _df.columns:
                    _low_col = _col.lower()
                    if _low_col == 'n' or ('sample' in _low_col and 'size' in _low_col):
                        _detected_n = int(_df[_col].max())
                        print(f"  N detected from column '{_col}': {_detected_n}")
                        break
                if _detected_n is None:
                    _file_id = os.path.basename(_harmonized_file).split('-')[0]
                    _gwas_dir = os.path.join("data", "gwas")
                    if os.path.exists(_gwas_dir):
                        for _f in os.listdir(_gwas_dir):
                            if _file_id in _f and _f.endswith(".yaml"):
                                with open(os.path.join(_gwas_dir, _f), 'r') as _fh:
                                    _detected_n = _find_n_deep(_yaml.safe_load(_fh))
                                    if _detected_n:
                                        print(f"  N detected from YAML ({_f}): {_detected_n}")
                                        break
                if _detected_n is None:
                    print("  No N column found — using study sample size N=807553")
                    _detected_n = 807553
                _df["N"] = _detected_n

            _hm3 = pd.read_csv(W_HM3_SNPLIST, sep="\t")[["SNP", "A1", "A2"]]
            _df  = _df.merge(_hm3, on="SNP", suffixes=("", "_hm3"))
            _df  = _df[
                (_df["A1"] == _df["A1_hm3"]) | (_df["A1"] == _df["A2_hm3"])
            ].drop(columns=["A1_hm3", "A2_hm3"])

            _strand_ambig = _df.apply(
                lambda r: set([r["A1"], r["A2"]]) in [{"A","T"}, {"C","G"}], axis=1
            )
            _df = _df[~_strand_ambig]

            _keep = [c for c in ["SNP","A1","A2","Z","N"] if c in _df.columns]
            _out  = _df[_keep].dropna(subset=["SNP","A1","A2","Z"])
            _out.to_csv(SUMSTATS_FILE, sep="\t", index=False, compression="gzip")
            print(f"  Written {len(_out):,} SNPs to {SUMSTATS_FILE}")
            print(f"  Mean |Z|: {_out['Z'].abs().mean():.3f}")
    return


@app.cell
def _(mo):
    mo.md("""
    ## 7. Generate cell-type-specific binary annotations (BED -> .annot.gz)
    """)
    return


@app.cell
def _(
    S3_BASE,
    all_cell_types,
    cell_type_beds,
    ldsc27_path,
    os,
    python27_path,
    subprocess,
):
    os.makedirs("data/annotations", exist_ok=True)
    _env = os.environ.copy()
    _env["PATH"] = f"{ldsc27_path}/bin:" + _env.get("PATH", "")

    for _ct in all_cell_types:
        _all_exist = all(
            os.path.exists(f"data/annotations/{_ct}.{_ch}.annot.gz")
            for _ch in range(1, 23)
        )
        if _all_exist:
            print(f"  {_ct}: all annotations exist, skipping")
            continue

        print(f"\nProcessing {_ct}...")
        _bed = cell_type_beds[_ct]
        for _ch in range(1, 23):
            _out = f"data/annotations/{_ct}.{_ch}.annot.gz"
            if os.path.exists(_out):
                continue
            _r = subprocess.run(
                [
                    python27_path, "tools/ldsc/make_annot.py",
                    "--bed-file",   _bed,
                    "--bimfile",    f"{S3_BASE}/data/reference/GRCh38/plink_files/1000G.EUR.hg38.{_ch}.bim",
                    "--annot-file", _out,
                ],
                capture_output=True, text=True, env=_env,
            )
            if _r.returncode != 0:
                print(f"  ERROR chr {_ch}: {_r.stderr}")
                raise subprocess.CalledProcessError(_r.returncode, _r.args)
            print(f"  chr{_ch}", end=" ", flush=True)
        print(f"\n  {_ct} done")

    print("\nAll annotations generated")
    return


@app.cell
def _(mo):
    mo.md("""
    ## 8. Calculate LD scores (HapMap3 SNPs only)
    """)
    return


@app.cell
def _(
    HM3_NO_MHC_LIST,
    S3_BASE,
    all_cell_types,
    concurrent,
    multiprocessing,
    os,
    python27_path,
    subprocess,
):
    os.makedirs("data/ldscores", exist_ok=True)

    def _calc_ld(args):
        ct, ch = args
        _dir  = f"data/ldscores/{ct}"
        os.makedirs(_dir, exist_ok=True)
        _out  = f"{_dir}/{ct}.{ch}.l2.ldscore.gz"
        if os.path.exists(_out):
            return f"[{ct}] chr{ch} exists"
        try:
            subprocess.run(
                [
                    python27_path, "tools/ldsc/ldsc.py",
                    "--l2",
                    "--bfile",      f"{S3_BASE}/data/reference/GRCh38/plink_files/1000G.EUR.hg38.{ch}",
                    "--ld-wind-cm", "1.0",
                    "--annot",      f"data/annotations/{ct}.{ch}.annot.gz",
                    "--thin-annot",
                    "--print-snps", HM3_NO_MHC_LIST,
                    "--out",        f"{_dir}/{ct}.{ch}",
                ],
                check=True, capture_output=True,
            )
            return f"[{ct}] chr{ch} done"
        except subprocess.CalledProcessError as e:
            return f"ERROR [{ct}] chr{ch}: {e.stderr}"

    _tasks      = [(ct, ch) for ct in all_cell_types for ch in range(1, 23)]
    _max_workers = min(multiprocessing.cpu_count() - 1, 6)
    print(f"Running LD score calculation ({_max_workers} workers, {len(_tasks)} tasks)...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=_max_workers) as _ex:
        for _res in _ex.map(_calc_ld, _tasks):
            print(_res)

    print("\nAll LD scores calculated")
    return


@app.cell
def _(mo):
    mo.md("""
    ## 9. Create CTS reference file
    """)
    return


@app.cell
def _(CTS_FILE, all_cell_types, os):
    os.makedirs("new_results", exist_ok=True)

    COMPLETED_CELL_TYPES = []
    for _ct in all_cell_types:
        _complete = all(
            os.path.exists(f"data/ldscores/{_ct}/{_ct}.{_ch}.l2.ldscore.gz")
            for _ch in range(1, 23)
        )
        if _complete:
            COMPLETED_CELL_TYPES.append(_ct)
        else:
            _missing = [c for c in range(1, 23)
                        if not os.path.exists(f"data/ldscores/{_ct}/{_ct}.{c}.l2.ldscore.gz")]
            print(f"  {_ct}: INCOMPLETE — missing chr {_missing}")

    with open(CTS_FILE, "w") as _f:
        for _ct in COMPLETED_CELL_TYPES:
            if _ct == "all_merged_cCREs":
                continue
            _f.write(f"{_ct}\tdata/ldscores/{_ct}/{_ct}.\n")

    print(f"\nCTS file written: {CTS_FILE}  ({len(COMPLETED_CELL_TYPES)} cell types)")
    return


@app.cell
def _(mo):
    mo.md("""
    ## 10. Run LDSC cell-type-specific heritability analysis
    """)
    return


@app.cell
def _(
    CTS_FILE,
    RESULTS_PREFIX,
    S3_BASE,
    SUMSTATS_FILE,
    os,
    python27_path,
    subprocess,
):
    if not os.path.exists(SUMSTATS_FILE):
        print(f"Skipping — sumstats not found: {SUMSTATS_FILE}")
    elif not os.path.exists(CTS_FILE):
        print(f"Skipping — CTS file not found: {CTS_FILE}")
    else:
        os.makedirs("new_results", exist_ok=True)
        print("Running LDSC CTS analysis...")
        subprocess.run(
            [
                python27_path, "tools/ldsc/ldsc.py",
                "--h2-cts",         SUMSTATS_FILE,
                "--ref-ld-chr",     (
                    f"{S3_BASE}/data/OSF/baseline_v1.2/baseline.,"
                    "data/ldscores/all_merged_cCREs/all_merged_cCREs."
                ),
                "--ref-ld-chr-cts", CTS_FILE,
                "--w-ld-chr",       f"{S3_BASE}/data/OSF/weights/weights.hm3_noMHC.",
                "--out",            RESULTS_PREFIX,
            ],
            check=False,
        )
        print("LDSC CTS analysis complete")
    return


@app.cell
def _(ollama, strong_research_prompt):

    response = ollama.chat(
        model='gemma4-bio',
        messages=[
            {'role': 'user', 'content': strong_research_prompt},
            {'role': 'assistant', 'content': '<think></think>\n'} 
        ]
    )

    gemma_raw = response['message']['content']
    print(gemma_raw)
    return


@app.cell
def _(genai, re, strong_research_prompt, time):
    _model = genai.GenerativeModel('gemini-3.5-flash')
    _raw_output = None
    for _attempt in range(3):
        try:
            _response = _model.generate_content(strong_research_prompt)
            _raw_output = _response.text
            break 
        except Exception as _e:
            if "429" in str(_e):
                print(f"RATE LIMIT: Pausing 40s for server cooldown...")
                time.sleep(40)
            else:
                print(f"CRITICAL ERROR: {_e}")
                break
    gemini_estimates = []
    if _raw_output:
        _names = re.findall(r'"cell_type":\s*"([^"]+)"', _raw_output)
        _reasons = re.findall(r'"(?:technical_reasoning|reason)":\s*"([^"]+)"', _raw_output)

        for _n, _r in zip(_names, _reasons):
            gemini_estimates.append({
                "cell_type": _n,
                "reasoning": _r
            })
    gemini_final_list = gemini_estimates[:3]
    if not gemini_final_list:
        print("No valid cell types identified. Check raw output:")
        print(_raw_output)
    else:
        for _i, _item in enumerate(gemini_final_list, 1):
            print(f"RANK {_i}: {_item['cell_type']}")
            print(f"reason: {_item['reasoning']}\n")
    return


@app.cell
def _(RESULTS_PREFIX, os, pd):
    from statsmodels.stats.multitest import fdrcorrection as _fdr

    _results_file = f"{RESULTS_PREFIX}.cell_type_results.txt"

    if not os.path.exists(_results_file):
        print(f"Results file not found: {_results_file}")
        ranked = None
    else:
        _results = pd.read_csv(_results_file, sep="\t")
        _, _results["FDR"] = _fdr(_results["Coefficient_P_value"].fillna(1))
        ranked = _results.sort_values("Coefficient_P_value")

        ranked.to_csv(f"{RESULTS_PREFIX}_ranked.csv", index=False)
        ranked.to_csv(f"{RESULTS_PREFIX}_ranked.txt", sep="\t", index=False)

        print("\nTop enriched cell types:\n")
        print(
            ranked[["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value", "FDR"]]
            .head(15)
            .to_string(index=False)
        )
        print(f"\nFDR < 0.05: {(ranked['FDR'] < 0.05).sum()} cell types")
        print(f"Saved to  : {RESULTS_PREFIX}_ranked.csv")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    jhj
    """)
    return


@app.cell
def _(ThreadPoolExecutor, as_completed, os, requests, time):
    PROJECT_ID = "86upf"
    TARGET_PATH = ["LDSC_hg38", "summary_statistics", "AlkesGroup"]
    DOWNLOAD_DIR = "data/gwas"
    SPECIFIC_FILES = ["PASS_AtrialFibrillation_Nielsen2018.sumstats.gz"]

    def get_osf_files(url):
        items = []
        while url:
            response = requests.get(url).json()
            items.extend(response['data'])
            url = response['links'].get('next')
        return items

    def download_file(url, filename, retries=3):
        path = os.path.join(DOWNLOAD_DIR, filename)
        tmp_path = path + ".tmp"

        if os.path.exists(path):
            print(f"Skipping {filename}, already exists.")
            return filename, True

        for attempt in range(retries):
            try:
                response = requests.get(url, stream=True, timeout=60)
                response.raise_for_status()
                total_size = int(response.headers.get('content-length', 0))
                downloaded = 0

                with open(tmp_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=4 * 1024 * 1024):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)

                if total_size and downloaded != total_size:
                    raise ValueError(f"Size mismatch: got {downloaded}, expected {total_size}")

                os.rename(tmp_path, path)
                print(f"Downloaded {filename} ({downloaded / 1e6:.1f} MB)")
                return filename, True

            except Exception as e:
                print(f"Attempt {attempt+1}/{retries} failed for {filename}: {e}")
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)

        return filename, False

    def run_dataset_download():
        if not SPECIFIC_FILES:
            print("ERROR: SPECIFIC_FILES list is empty.")
            return

        os.makedirs(DOWNLOAD_DIR, exist_ok=True)

        print(f"Connecting to OSF Project: {PROJECT_ID}...")
        api_url = f"https://api.osf.io/v2/nodes/{PROJECT_ID}/files/osfstorage/"
        current_items = get_osf_files(api_url)

        for folder_name in TARGET_PATH:
            for item in current_items:
                if item['attributes']['kind'] == 'folder' and item['attributes']['name'] == folder_name:
                    api_url = item['relationships']['files']['links']['related']['href']
                    current_items = get_osf_files(api_url)
                    break

        url_map = {
            item['attributes']['name']: item['links']['download']
            for item in current_items
            if item['attributes']['kind'] == 'file' and item['attributes']['name'] in SPECIFIC_FILES
        }

        missing = set(SPECIFIC_FILES) - set(url_map.keys())
        if missing:
            print(f"WARNING: files not found on server: {missing}")

        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {
                executor.submit(download_file, url, name): name
                for name, url in url_map.items()
            }
            for future in as_completed(futures):
                name, success = future.result()
                if not success:
                    print(f"FAILED: {name}")

    run_dataset_download()
    return


if __name__ == "__main__":
    app.run()
