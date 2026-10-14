import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    return (mo,)


@app.cell
def __(mo):
    mo.md(r"""
    # Drosophila Female Lifespan Fine-Mapping
    Configure inputs and run a selected stage of the GWAS and fine-mapping pipeline. The original command-line pipeline remains available beside this interface.
    """)
    return


@app.cell
def __():
    from pathlib import Path
    import subprocess
    import sys

    project_root = Path(__file__).resolve().parents[1]
    pipeline_script = project_root / "notebooks" / "drosophila_female_lifespan_finemapping_pipeline.py"
    defaults = {
        "base_dir": project_root,
        "pheno": project_root / "data/gwas/lifespan_female.pheno",
        "glm_dir": project_root / "data/gwas/tmp",
        "qc_dir": project_root / "data/gwas/tmp/qc",
        "eigenvec": project_root / "data/gwas/tmp/dgrp_pca.eigenvec",
        "bfile": project_root / "data/gwas/tmp/merged_qc",
        "gtf": project_root / "data/genes/Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz",
        "enhancer_dir": project_root / "data/enhancers/dm",
        "out_dir": project_root / "data/finemap/female",
    }
    return defaults, pipeline_script, project_root, subprocess, sys


@app.cell
def __(mo):
    mo.md(r"""
    ## Inputs and Run Settings
    """)
    return


@app.cell
def __(defaults, mo):
    stage_input = mo.ui.dropdown(
        options=["all", "prepare", "plots", "gwas", "cojo", "susie", "genes", "enhancers"],
        value="all", label="Pipeline stage",
    )
    base_dir_input = mo.ui.text(value=str(defaults["base_dir"]), label="Project directory")
    pheno_input = mo.ui.text(value=str(defaults["pheno"]), label="Phenotype file")
    glm_dir_input = mo.ui.text(value=str(defaults["glm_dir"]), label="GWAS working directory")
    qc_dir_input = mo.ui.text(value=str(defaults["qc_dir"]), label="QC genotype directory")
    eigenvec_input = mo.ui.text(value=str(defaults["eigenvec"]), label="PCA covariates")
    bfile_input = mo.ui.text(value=str(defaults["bfile"]), label="Merged PLINK genotype prefix")
    gtf_input = mo.ui.text(value=str(defaults["gtf"]), label="Drosophila GTF annotation")
    enhancer_dir_input = mo.ui.text(value=str(defaults["enhancer_dir"]), label="Enhancer BED directory")
    out_dir_input = mo.ui.text(value=str(defaults["out_dir"]), label="Fine-mapping output directory")
    gcta_input = mo.ui.text(value="gcta64", label="GCTA executable")
    plink_input = mo.ui.text(value="plink", label="PLINK executable")
    plink2_input = mo.ui.text(value="plink2", label="PLINK2 executable")
    force_input = mo.ui.checkbox(value=False, label="Regenerate cached outputs")
    mo.vstack([
        stage_input, base_dir_input, pheno_input, glm_dir_input, qc_dir_input, eigenvec_input,
        bfile_input, gtf_input, enhancer_dir_input, out_dir_input,
        gcta_input, plink_input, plink2_input, force_input,
    ])
    return (base_dir_input, bfile_input, eigenvec_input, enhancer_dir_input, force_input,
            gcta_input, gtf_input, glm_dir_input, out_dir_input, pheno_input, plink2_input,
            plink_input, qc_dir_input, stage_input)


@app.cell
def __(mo):
    run_button = mo.ui.run_button(label="Run selected fine-mapping stage")
    run_button
    return (run_button,)


@app.cell
def __(base_dir_input, bfile_input, eigenvec_input, enhancer_dir_input, force_input,
      gcta_input, gtf_input, glm_dir_input, mo, out_dir_input, pheno_input, pipeline_script,
      plink2_input, plink_input, qc_dir_input, run_button, stage_input, subprocess, sys):
    mo.stop(not run_button.value, mo.md("Set the inputs above, choose a stage, then press the run button."))
    command = [
        sys.executable, str(pipeline_script), "--stage", stage_input.value,
        "--base-dir", base_dir_input.value,
        "--pheno", pheno_input.value,
        "--glm-dir", glm_dir_input.value,
        "--qc-dir", qc_dir_input.value,
        "--eigenvec", eigenvec_input.value,
        "--bfile", bfile_input.value,
        "--gtf", gtf_input.value,
        "--enhancer-dir", enhancer_dir_input.value,
        "--out-dir", out_dir_input.value,
        "--gcta-bin", gcta_input.value,
        "--plink-bin", plink_input.value,
        "--plink2-bin", plink2_input.value,
    ]
    if force_input.value:
        command.append("--force")
    print("Running:", " ".join(command), flush=True)
    try:
        subprocess.run(command, cwd=pipeline_script.parents[1], check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Fine-mapping stage '{stage_input.value}' failed with exit code {error.returncode}.") from error
    print(f"Fine-mapping stage '{stage_input.value}' completed successfully.")
    return


if __name__ == "__main__":
    app.run()
