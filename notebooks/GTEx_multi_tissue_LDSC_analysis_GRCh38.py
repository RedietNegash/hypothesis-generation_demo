import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():

    import marimo as mo
    import urllib.request
    import os
    import subprocess
    import pandas as pd
    import tarfile
    from pathlib import Path
    import glob
    import json
    return mo, urllib, os, subprocess, pd, tarfile, Path, glob, json


@app.cell
def __(mo):
    mo.md("""
# LDSC Multi-Tissue Analysis – GTEx Gene Sets

This notebook performs tissue-specific heritability analysis using LDSC (LD Score Regression)  
with GTEx tissue-specific gene sets. The analysis identifies which tissues contribute most  
to the heritability of complex traits.

All analyses are performed using **GRCh38 / hg38** coordinates.

## Analysis Overview:
1. Process tissue-specific top gene lists
2. Create genomic annotations for each tissue
3. Calculate tissue-specific LD scores
4. Run partitioned heritability analysis
""")
    return


@app.cell
def __(mo):
    mo.md("""
## 0. Setup: Download and configure LDSC

This cell will:
1. Create a Python 2.7 conda environment for LDSC
2. Install LDSC and its dependencies
3. Install BEDTools (required for annotation generation)
""")
    return


@app.cell
def __(Path, subprocess, os, json):
    TOOLS_DIR = Path("tools")
    LDSC_DIR = TOOLS_DIR / "ldsc"
    TOOLS_DIR.mkdir(exist_ok=True)
    
 
    env_check = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True,
        text=True
    )
    ldsc_env_exists = "ldsc27" in env_check.stdout
    
    if not ldsc_env_exists:
        print("Creating Python 2.7 conda environment for LDSC...")
        subprocess.run([
            "conda", "create", "-n", "ldsc27", 
            "python=2.7", "-y"
        ], check=True)
    
  
    conda_prefix = subprocess.run(
        ["conda", "env", "list", "--json"],
        capture_output=True,
        text=True,
        check=True
    )
    envs = json.loads(conda_prefix.stdout)["envs"]
    ldsc27_path = [e for e in envs if "ldsc27" in e][0]
    
   
    bedtools_check = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "bedtools"), "--version"],
        capture_output=True,
        text=True
    )
    
    if bedtools_check.returncode != 0:
        print("Installing BEDTools in ldsc27 environment...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27",
            "-c", "bioconda", "bedtools", "-y"
        ], check=True)
    else:
        print(f" BEDTools already installed: {bedtools_check.stdout.strip()}")
    
 
    if not LDSC_DIR.exists():
        print("Cloning LDSC repository...")
        subprocess.run([
            "git", "clone",
            "https://github.com/bulik/ldsc.git",
            str(LDSC_DIR)
        ], check=True)
    
    check_numpy = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "python"), "-c", "import numpy"],
        capture_output=True
    )
    
    if check_numpy.returncode != 0:
        print("Installing LDSC dependencies in ldsc27 environment...")
        print("Installing OpenSSL 1.0...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "openssl=1.0.2", "-c", "conda-forge"
        ], check=True)
        
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "numpy", "scipy", "pandas", "bitarray", "-c", "conda-forge"
        ], check=True)
        
        print("Installing pybedtools and pysam...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "pybedtools", "pysam=0.15.3", "-c", "bioconda", "-c", "conda-forge"
        ], check=True)
        
        print("Dependencies installed!")
    
    print("\nVerifying BEDTools is accessible to pybedtools...")
    pybedtools_check = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "python"), "-c", 
         "from pybedtools import BedTool; import pybedtools.helpers as helpers; print('BEDTools path:', helpers.get_bedtools_path())"],
        capture_output=True,
        text=True
    )
    
    if pybedtools_check.returncode != 0:
        print("WARNING: pybedtools can't find BEDTools!")
        print("Error:", pybedtools_check.stderr)
    else:
        print("pybedtools can access BEDTools")
        print(pybedtools_check.stdout)
    
    ldsc_script = LDSC_DIR / "ldsc.py"
    subprocess.run(["chmod", "+x", str(ldsc_script)], check=True)
    
    python27_path = os.path.join(ldsc27_path, "bin", "python")
    
    print(f"\n LDSC environment ready!")
    print(f"Python 2.7 path: {python27_path}")
    print(f"LDSC path: {ldsc_script}")
    
    return (ldsc_script, python27_path, ldsc27_path)

app.cell
def __(mo):
    mo.md("## 1. Download GTEx gene sets, reference panels, and GWAS summary statistics")
    return


@app.cell
def __(os, urllib):

    os.makedirs("data/gtex/tissue_gene_sets", exist_ok=True)
    os.makedirs("data/reference", exist_ok=True)
    os.makedirs("data/gwas", exist_ok=True)
    os.makedirs("data/gtex/annot_files", exist_ok=True)
    os.makedirs("data/gtex/ldscores", exist_ok=True)
    
    print(" Data directories created")
    return


@app.cell
def __(mo):

    mo.md("""
### GTEx Tissue Gene Sets

**Important**: Place your GTEx tissue-specific gene lists in `data/gtex/tissue_gene_sets/`

Expected file format: `{tissue_name}_top10.txt` (or your chosen suffix)

Each file should contain:
- One gene per line (gene symbols or Ensembl IDs)
- Top genes expressed in that tissue
- Example tissues: Brain_Cortex, Liver, Heart, etc.

You can also modify the code below to specify custom file path.
""")
    return


@app.cell
def __(os, urllib):

    if not os.path.exists("data/reference/GRCh38.tgz"):
        print("\nDownloading GRCh38 reference with baseline LD scores (this may take several minutes)...")
        urllib.request.urlretrieve(
            "https://zenodo.org/records/10515792/files/GRCh38.tgz?download=1",
            "data/reference/GRCh38.tgz"
        )
        print(" GRCh38 reference downloaded")
    else:
        print("\n GRCh38 reference already downloaded")
    
    if not os.path.exists("data/gwas/AD_bellenguez_2022_hg38.tsv.gz"):
        print("Downloading example GWAS summary statistics (Alzheimer's disease)...")
        print("You can replace this with your trait of interest.")
        urllib.request.urlretrieve(
            "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz",
            "data/gwas/AD_bellenguez_2022_hg38.tsv.gz"
        )
        print(" GWAS data downloaded")
    else:
        print(" GWAS data already downloaded")
    
    print("\n All downloads complete!")
    return


@app.cell
def __(mo):

    mo.md("## 2. Extract reference LD panels and baseline LD scores")
    return


@app.cell  
def __(tarfile, os):
    print("Checking GRCh38.tgz contents...")
    
    if os.path.exists("data/reference/GRCh38.tgz"):
        with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
            members = _tar.getmembers()
            print(f"Archive contains {len(members)} items")
            print("\nTop-level structure:")
            seen = set()
            for m in members[:50]:  
                parts = m.name.split('/')
                if len(parts) > 1:
                    top = parts[0] + "/" + parts[1]
                    if top not in seen:
                        print(f"  {top}")
                        seen.add(top)
    return


@app.cell
def __(tarfile, os):

    print("Extracting GRCh38 reference files...")
    
    if not os.path.exists("data/reference/GRCh38"):
        try:
            print("  Verifying GRCh38.tgz integrity...")
            with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
                _tar.getmembers()
            print("  File integrity verified")
        except (EOFError, tarfile.ReadError) as e:
            print(f"\n  Error: GRCh38.tgz is corrupted!")
            print(f"  Please delete it and re-run: rm data/reference/GRCh38.tgz")
            raise
        
        print("  Extracting (this may take a few minutes)...")
        with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
            _tar.extractall("data/reference")
        print("GRCh38 reference extracted successfully")
    else:
        print("GRCh38 directory exists")
    
  
    nested_files = [
        ("data/reference/GRCh38/baselineLD_v2.2.tgz", "data/reference", "baselineLD_v2.2"),
        ("data/reference/GRCh38/plink_files.tgz", "data/reference/GRCh38", "1000G.EUR.hg38.1.bed"),
        ("data/reference/GRCh38/weights.tgz", "data/reference/GRCh38", "weights")
    ]
    
    for tar_file, extract_to, check_file in nested_files:
        check_path = os.path.join(extract_to, check_file)
        if os.path.exists(check_path):
            print(f"  {os.path.basename(tar_file)} already extracted")
            continue
            
        if os.path.exists(tar_file):
            print(f"  Extracting {os.path.basename(tar_file)}...")
            with tarfile.open(tar_file, "r:gz") as _tar:
                _tar.extractall(extract_to)
            print(f" {os.path.basename(tar_file)} extracted")
            
            if not os.path.exists(check_path):
                print(f" Warning: Expected file {check_path} not found after extraction")
        else:
            print(f"  Warning: {tar_file} not found")
    
    critical_file = "data/reference/GRCh38/1000G.EUR.hg38.1.bim"
    if os.path.exists(critical_file):
        print(f"\n All reference files ready! Verified: {critical_file}")
    else:
        print(f"\n ERROR: Critical file missing: {critical_file}")
        print("Checking what files exist in data/reference/GRCh38/:")
        if os.path.exists("data/reference/GRCh38"):
            files = os.listdir("data/reference/GRCh38")
            print(f"  Found {len(files)} files/directories")
            for _f in sorted(files)[:10]:  
                print(f"    - {_f}")
        else:
            print("  Directory doesn't exist!")
    
    return

@app.cell
def __(mo):
    """
    COMMIT: Add GWAS munging section header
    """
    mo.md("## 3. Munge GWAS summary statistics")
    return


@app.cell
def __(subprocess, os, python27_path):
    os.makedirs("data/munged", exist_ok=True)

    munged_file = "data/munged/AD_bellenguez_2022_hg38_munged.sumstats.gz"
    
    if os.path.exists(munged_file):
        print(" Munged GWAS file already exists, skipping munging step")
    else:
        print("\n" + "="*60)
        print("STEP 3: Munging GWAS summary statistics")
        print("="*60)
        print("\nNote: Modify column names below to match your GWAS format!")
        
        subprocess.run([
            python27_path, "tools/ldsc/munge_sumstats.py",
            "--sumstats", "data/gwas/AD_bellenguez_2022_hg38.tsv.gz",
            "--out", "data/munged/AD_bellenguez_2022_hg38_munged",
            "--a1", "effect_allele",
            "--a2", "other_allele",
            "--p", "p_value",
            "--snp", "variant_id",
            "--N-col", "n_total"
        ], check=True)
        
        print("\n GWAS munging complete")
    return

@app.cell
def __(mo):
    """
    COMMIT: Add tissue discovery section header
    """
    mo.md("""
## 4. Discover GTEx tissue gene sets

This step will automatically find all `*_top10.txt` files in the tissue gene sets directory.
You can modify the pattern to match your file naming convention.
""")
    return


@app.cell
def __(glob, os):
    gene_set_pattern = "data/gtex/tissue_gene_sets/*_top10.txt"
    gene_set_files = glob.glob(gene_set_pattern)
    
    if not gene_set_files:
        print(f"  No gene set files found matching: {gene_set_pattern}")
        print("\nPlease ensure your tissue gene lists are in the correct location:")
        print("  data/gtex/tissue_gene_sets/{tissue_name}_top10.txt")
        tissues = []
    else:
        tissues = []
        for filepath in gene_set_files:
            basename = os.path.basename(filepath)
            tissue_name = basename.replace("_top10.txt", "")
            tissues.append(tissue_name)
        
        print(f" Found {len(tissues)} tissue gene sets:")
        for tissue in sorted(tissues):
            filepath = f"data/gtex/tissue_gene_sets/{tissue}_top10.txt"
            try:
                with open(filepath) as f:
                    gene_count = len(f.readlines())
                print(f"  - {tissue} ({gene_count} genes)")
            except:
                print(f"  - {tissue}")
    
    return (tissues, gene_set_files)
@app.cell
def __(mo):
    mo.md("""
## 5. Generate tissue-specific genomic annotations

For each tissue and chromosome, this creates binary annotations indicating which SNPs  
fall within genes that are highly expressed in that tissue.

""")
    return

@app.cell
def __(mo):
    mo.md("""
### Gene Coordinate Mapping Strategy

To create annotations

**Option 1**: Use a pre-built gene annotation file (GTF/GFF from Ensembl/GENCODE)
**Option 2**: Create BED files manually for each tissue's gene list
**Option 3**: Use biomaRt or similar tool to fetch coordinates programmatically

For this notebook, we'll demonstrate **Option 3** using a simple approach.
You may need to install additional packages or provide your own coordinate files.
""")
    return


@app.cell
def __(pd, subprocess, os, tissues, python27_path, ldsc27_path):

    print("\n" + "="*60)
    print("STEP 5: Generating tissue-specific annotations")
    print("="*60)
    
    env = os.environ.copy()
    env["PATH"] = f"{ldsc27_path}/bin:" + env.get("PATH", "")
    
    
    print("\n IMPORTANT: Gene coordinate mapping required!")
    print("This step requires gene coordinates. You need to:")
    print("1. Obtain a gene annotation file (GTF/GFF)")
    print("2. Map gene symbols to genomic coordinates")
    print("3. Create BED files for each tissue")
    print("\nSkipping annotation generation for now.")
    print("Please implement gene coordinate mapping based on your data source.")
    
    return
@app.cell
def __(mo):

    mo.md("""
## 6. Calculate LD scores for each tissue and chromosome

This step computes LD scores using the tissue-specific annotations.
This is the computational bottleneck and may take 10-30 minutes per tissue.
""")
    return


@app.cell
def __(subprocess, os, tissues, python27_path):
    os.makedirs("data/gtex/ldscores", exist_ok=True)

    print("\n" + "="*60)
    print("STEP 6: Calculating LD scores")
    print("This matches your original script:")
    print("  --l2 --ld-wind-cm 1 --thin-annot")
    print("="*60)
    
    for tissue in tissues:
        print(f"\nProcessing {tissue}...")
        os.makedirs(f"data/gtex/ldscores/{tissue}", exist_ok=True)
        _all_exist = all(
            os.path.exists(f"data/gtex/ldscores/{tissue}/{tissue}.{chrom}.l2.ldscore.gz")
            for chrom in range(1, 23)
        )
        
        if _all_exist:
            print(f" All {tissue} LD scores already exist, skipping")
            continue

        for chrom in range(1, 23):
            ldscore_file = f"data/gtex/ldscores/{tissue}/{tissue}.{chrom}.l2.ldscore.gz"
            
            if os.path.exists(ldscore_file):
                print(f"  Chromosome {chrom}  (exists)", end=" ", flush=True)
                continue

            annot_file = f"data/gtex/annot_files/{tissue}.{chrom}.annot.gz"
            if not os.path.exists(annot_file):
                print(f"  Chromosome {chrom}  (annotation missing)")
                continue
                
            print(f"  Chromosome {chrom}...", end=" ", flush=True)
            subprocess.run([
                python27_path, "tools/ldsc/ldsc.py",
                "--l2",
                "--bfile", f"data/reference/GRCh38/plink_files/1000G.EUR.hg38.{chrom}",
                "--ld-wind-cm", "1",
                "--annot", annot_file,
                "--thin-annot",
                "--out", f"data/gtex/ldscores/{tissue}/{tissue}.{chrom}"
            ], check=True, capture_output=True)
            
            print("✓")
        
        print(f"  {tissue} complete")
    
    print("\n All LD scores calculated")
    return
@app.cell
def __(mo):
    mo.md("## 7. Create tissue-specific CTS reference file")
    return


@app.cell
def __(os, tissues):
    os.makedirs("results", exist_ok=True)
    
    cts_file = "data/gtex_tissues.cts"
    with open(cts_file, "w") as f:
        for tissue in tissues:
            f.write(f"{tissue}    data/gtex/ldscores/{tissue}/{tissue}.\n")
    
    print(f"CTS reference file created: {cts_file}")
    print(f"Contains {len(tissues)} tissues")

    print("\nContents:")
    with open(cts_file) as f:
        for line in f:
            print(f"  {line.strip()}")
    
    return (cts_file,)
@app.cell
def __(mo):
    mo.md("""
## 8. Run LDSC tissue-specific heritability analysis

This performs partitioned heritability analysis to identify which tissues  
contribute significantly to trait heritability.
""")
    return


@app.cell
def __(subprocess, python27_path, cts_file):
    print("\n" + "="*60)
    print("STEP 8: Running tissue-specific heritability analysis")
    print("="*60 + "\n")
    
    output_prefix = "results/GTEx_TissueSpecific"
    
    if os.path.exists(f"{output_prefix}.cell_type_results.txt"):
        print("Results file already exists, skipping analysis")
    else:
        subprocess.run([
            python27_path, "tools/ldsc/ldsc.py",
            "--h2-cts", "data/munged/AD_bellenguez_2022_hg38_munged.sumstats.gz",
            "--ref-ld-chr", "data/reference/baselineLD_v2.2/baselineLD.",
            "--ref-ld-chr-cts", cts_file,
            "--w-ld-chr", "data/reference/GRCh38/weights/weights.hm3_noMHC.",
            "--out", output_prefix
        ], check=True)
        
        print("\n Tissue-specific analysis complete")
    
    return (output_prefix,)

@app.cell
def __(mo):
    mo.md("## 9. Rank tissues by heritability enrichment significance")
    return


@app.cell
def __(pd, output_prefix):
    print("\n" + "="*60)
    print("STEP 9: Ranking tissues by significance")
    print("="*60 + "\n")
    
    results_file = f"{output_prefix}.cell_type_results.txt"
    
    if not os.path.exists(results_file):
        print(f" Results file not found: {results_file}")
        print("Please run the heritability analysis first (Step 8)")
    else:
        results = pd.read_csv(results_file, sep="\t")
        ranked = results.sort_values("Coefficient_P_value")
        
        output_csv = "results/GTEx_TissueSpecific_ranked.csv"
        ranked.to_csv(output_csv, index=False)
        
        print("Tissues ranked by heritability enrichment p-value:\n")
        print(ranked[["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"]].to_string(index=False))
        
        print(f"\n✓ Results saved to {output_csv}")
        significant = ranked[ranked["Coefficient_P_value"] < 0.05]
        if len(significant) > 0:
            print(f"\n Found {len(significant)} tissue(s) with significant enichment (p < 0.05):")
            for idx, row in significant.iterrows():
                print(f"  - {row['Name']}: p = {row['Coefficient_P_value']:.2e}")
        
        ranked_results = ranked
    
    return (ranked, ranked_results)


@app.cell
def __(mo):
    mo.md("""
## Next Steps

### Interpreting Results
- **Coefficient**: Heritability enrichment in each tissue
- **Coefficient_P_value**: Statistical significance of enrichment
- **Lower p-values** indicate stronger tissue-specific contribution to trait heritability


### Gene Coordinate Mapping
The annotation generation step (Step 5) requires gene coordinates. You can:
- Download GENCODE annotations: https://www.gencodegenes.org/
- Use Ensembl BioMart: https://www.ensembl.org/biomart/
- Provide your own BED files with gene coordinates

### References
- LDSC: Finucane et al. (2015) Nature Genetics
- Baseline LD Model: Gazal et al. (2017) Nature Genetics
- GTEx: GTEx Consortium (2020) Science
""")
    return


if __name__ == "__main__":
    app.run() 