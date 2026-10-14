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
    # Drosophila Female Lifespan MAGMA
    Configure inputs and run a selected stage of the MAGMA gene and gene-set analysis. The original command-line pipeline remains available beside this interface.
    """)
    return


@app.cell
def __():
    from pathlib import Path
    import sys
    import subprocess

    project_root = Path(__file__).resolve().parents[1]
    pipeline_script = project_root / "notebooks" / "drosophila_female_lifespan_magma_pipeline.py"
    defaults = {
        "sumstats": project_root / "data/gwas/tmp/lifespan_female_gwas_sumstats.tsv",
        "gtf": project_root / "data/genes/Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz",
        "bfile": project_root / "data/gwas/tmp/merged_qc",
        "magma_bin": project_root / "tools/magma/magma",
        "out_dir": project_root / "results/magma_female",
        "go_gaf": project_root / "results/magma_female/go_inputs/DROME-mod.gaf.gz",
        "go_obo": project_root / "results/magma_female/go_inputs/go-basic.obo",
        "gene_sets": project_root / "results/magma_female/go_biological_process.sets.txt",
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
    sumstats_input = mo.ui.text(value=str(defaults["sumstats"]), label="GWAS summary statistics")
    gtf_input = mo.ui.text(value=str(defaults["gtf"]), label="Drosophila GTF annotation")
    bfile_input = mo.ui.text(value=str(defaults["bfile"]), label="PLINK genotype prefix")
    magma_bin_input = mo.ui.text(value=str(defaults["magma_bin"]), label="MAGMA executable")
    out_dir_input = mo.ui.text(value=str(defaults["out_dir"]), label="Output directory")
    stage_input = mo.ui.dropdown(
        options=["all", "inputs", "install", "validate", "prepare", "annotate", "genes", "sets", "report"],
        value="all", label="Pipeline stage",
    )
    gene_source_input = mo.ui.dropdown(
        options=["No gene-set analysis", "Download Gene Ontology", "Use gene-set file", "Use GAF and OBO files"],
        value="Download Gene Ontology", label="Gene-set source",
    )
    gene_sets_input = mo.ui.text(value=str(defaults["gene_sets"]), label="Gene-set file")
    gaf_input = mo.ui.text(value=str(defaults["go_gaf"]), label="Gene Ontology GAF file")
    obo_input = mo.ui.text(value=str(defaults["go_obo"]), label="Gene Ontology OBO file")
    force_input = mo.ui.checkbox(value=False, label="Rebuild cached outputs")
    mo.vstack([
        sumstats_input, gtf_input, bfile_input, magma_bin_input, out_dir_input,
        stage_input, force_input, gene_source_input, gene_sets_input, gaf_input, obo_input,
    ])
    return (bfile_input, force_input, gaf_input, gene_sets_input, gene_source_input,
            magma_bin_input, obo_input, out_dir_input, stage_input, sumstats_input, gtf_input)


@app.cell
def __(mo):
    run_button = mo.ui.run_button(label="Run selected MAGMA stage")
    run_button
    return (run_button,)


@app.cell
def __(bfile_input, defaults, force_input, gaf_input, gene_sets_input, gene_source_input,
      magma_bin_input, mo, obo_input, out_dir_input, pipeline_script, run_button,
      stage_input, subprocess, sumstats_input, gtf_input, sys):
    mo.stop(not run_button.value, mo.md("Set the inputs above, choose a stage, then press the run button."))
    command = [
        sys.executable, str(pipeline_script),
        "--sumstats", sumstats_input.value,
        "--gtf", gtf_input.value,
        "--bfile", bfile_input.value,
        "--magma-bin", magma_bin_input.value,
        "--out-dir", out_dir_input.value,
        "--stage", stage_input.value,
    ]
    if force_input.value:
        command.append("--force")
    if gene_source_input.value == "Download Gene Ontology":
        command.append("--download-go")
    elif gene_source_input.value == "Use gene-set file":
        command.extend(["--gene-sets", gene_sets_input.value])
    elif gene_source_input.value == "Use GAF and OBO files":
        command.extend(["--go-gaf", gaf_input.value, "--go-obo", obo_input.value])
    print("Running:", " ".join(command), flush=True)
    try:
        subprocess.run(command, cwd=pipeline_script.parents[1], check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"MAGMA stage '{stage_input.value}' failed with exit code {error.returncode}.") from error
    print(f"MAGMA stage '{stage_input.value}' completed successfully.")
    return


if __name__ == "__main__":
    app.run()
