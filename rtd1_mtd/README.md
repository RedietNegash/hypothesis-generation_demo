# RTD-1 / MTD transcriptomics

Branch: `rtd1_mtd_transcriptomics`. Scope: **Specific Objectives 1 and 2 only**—collect/prepare biological datasets, analyse expression, identify cell populations provisionally, and examine pathway associations.

Open [the Jupyter notebook](RTD1_MTD_objectives_1_2.ipynb) for a guided walkthrough with saved results, plots, source checks, and optional full rerun cells. Its default mode reviews existing DESeq2/Scanpy fits and recomputes summary displays and validation.

For notebook execution, install `requirements-notebook.txt` into the project environment and select that environment as the notebook kernel. This adds Jupyter kernel and execution dependencies to the analysis requirements.

Start with the [analysis and findings](results/analysis_summary.md), [dataset inventory](data/dataset_inventory.tsv), and [validation results](results/validation.json). Existing repository projects are outside this directory and were not changed. No network integration, peptide-target ranking, AI mechanism prediction, or RTD-1/MTD comparison is implemented.

## Inputs and interpretation

- The [2019 rat arthritis supplement](https://doi.org/10.6084/m9.figshare.8855588) supplies processed gene lists, not a count matrix. Original values and disagreements with paper-level totals/directions are preserved. Enrichment is conditional on the published disease-gene list; it is not a reproduction of whole-transcriptome enrichment.
- [GSE278567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278567) supplies human THP-1 bulk estimated counts. Six DESeq2 dose-versus-vehicle contrasts are run separately within basal/LPS conditions. Fractional E/M counts are rounded explicitly, which makes this an exploratory reanalysis.
- [GSE269120](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269120) supplies rat control/clodronate synovial single-cell matrices. Scanpy produces QC, clusters, marker evidence, and descriptive pathways. These are not peptide-treated samples. One library per condition precludes replicate-based condition testing. Cell labels are provisional, with mixed/unresolved clusters retained.
- Cytokine/protein/MTD material is inventoried with availability gaps. No quantitative measurements were invented from figures and no author contact was made.

## Reproduce

Run from the repository root. Python 3.10 and system R 4.6 were used. Full R versions are in [sessionInfo](results/bulk/sessionInfo.txt); Python dependencies are pinned in `requirements.lock.txt`. R installation is isolated to `.Rlib`, and data/cache/environment files are ignored by Git.

```bash
uv venv --python 3.10 rtd1_mtd/.venv
uv pip install --python rtd1_mtd/.venv/bin/python -r rtd1_mtd/requirements.lock.txt
/usr/bin/Rscript rtd1_mtd/scripts/install_r.R
rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/download.py
rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/prepare.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 /usr/bin/Rscript rtd1_mtd/scripts/bulk_deseq2.R
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 /usr/bin/Rscript rtd1_mtd/scripts/bulk_deseq2.R --block-sensitivity
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/singlecell.py
rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/pathways.py
rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/report.py
rtd1_mtd/.venv/bin/python rtd1_mtd/scripts/validate.py
```

The download script reuses existing files and downloads missing files atomically. The committed source manifest records SHA256 values of this run, including Reactome release 98. Reactome `current` URLs can change: for an exact rerun, use cached files matching the manifest. Downloading a new release requires regenerating results and reviewing mapping changes. `download.py --manifest-only` records local source checksums without network access; `validate.py` verifies files against the recorded manifest.

The single-cell count object is stored locally at `data/processed/rat_synovium.h5ad` and is excluded from Git. Cluster labels in `data/cluster_annotation.tsv` are manually curated for the recorded clustering; they must be reviewed if data, settings, or dependency versions change. The R installer pins the Bioconductor repository release but resolves CRAN dependencies at installation time; the recorded `sessionInfo.txt` is the authoritative environment record for this run.

## Outputs

- `data/processed/`: prepared rat gene lists, human counts/sample metadata/QC, single-cell object, and species-specific Reactome mappings.
- `results/bulk/`: six DESeq2 result tables, contrast summary, sample QC/PCA, methods, and R session details.
- `results/bulk_sensitivity/`: six candidate-block contrasts and comparison with the primary model; replicate labels are not verified experimental blocks.
- `results/singlecell/`: QC, cluster markers, marker panels, cell annotations, descriptive cell counts, UMAP and dot plots. Cluster marker p values are cell-level descriptive statistics, not treatment tests.
- `results/pathways/`: per-list enrichment tables with all tested pathways, annotation coverage audit, and methods. No enrichment implies direct peptide binding or causal pathway activation.
- `results/analysis_summary.md`: findings, figures, source discrepancies, and limitations.
- `results/validation.json`: source integrity and scientific input/output consistency checks.
