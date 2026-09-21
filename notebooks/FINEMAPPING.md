# Female lifespan fine-mapping pipeline

The workflow has three stages:

1. `drosophila_female_lifespan_finemapping_pipeline.py --stage cojo` selects independent GWAS signals with GCTA-COJO and extracts their loci.
2. `drosophila_female_lifespan_finemapping_pipeline.py --stage susie` fine-maps each locus with SuSiE-RSS.
3. `fly_finemapping_gene_mapping_female.py` maps selected variants to genes.

Stages 1 and 2 are integrated in the unified pipeline.

## Stage 1 inputs

Place the chromosome-level PLINK 2 association results in `data/gwas/tmp`:

```text
lifespan_2L.S18_1537_F.glm.linear
lifespan_2R.S18_1537_F.glm.linear
lifespan_3L.S18_1537_F.glm.linear
lifespan_3R.S18_1537_F.glm.linear
lifespan_4.S18_1537_F.glm.linear
lifespan_X.S18_1537_F.glm.linear
```

Each association file must contain these columns:

```text
#CHROM POS ID A1 OMITTED A1_FREQ BETA SE P OBS_CT
```

The same directory must contain the QC-filtered PLINK LD reference:

```text
merged_qc.bed
merged_qc.bim
merged_qc.fam
```

The BIM chromosome field may use `2L`, `2R`, `3L`, `3R`, `4`, and `X`. Stage 1 remaps these labels to numeric GCTA chromosomes.

## Software

The Python environment must provide NumPy, pandas, and rpy2. R must provide the `susieR` package. GCTA and PLINK must either be available on `PATH` or supplied with `--gcta-bin` and `--plink-bin`.

Record the software versions used for a run:

```bash
python --version
gcta64 --version
plink --version
R --quiet --no-save -e 'packageVersion("susieR")'
```

## Run the pipeline

Run both integrated stages from the repository root:

```bash
python notebooks/drosophila_female_lifespan_finemapping_pipeline.py --stage all
```

Run a single stage:

```bash
python notebooks/drosophila_female_lifespan_finemapping_pipeline.py \
  --stage cojo \
  --gcta-bin /absolute/path/to/gcta64

python notebooks/drosophila_female_lifespan_finemapping_pipeline.py \
  --stage susie \
  --plink-bin /absolute/path/to/plink
```

Existing outputs are reused by default. Use `--force` after changing inputs or analysis settings:

```bash
python notebooks/drosophila_female_lifespan_finemapping_pipeline.py \
  --stage all \
  --gcta-bin /absolute/path/to/gcta64 \
  --plink-bin /absolute/path/to/plink \
  --force
```

## Analysis defaults

| Setting | Value |
| --- | ---: |
| Phenotype | `S18_1537_F` |
| Minimum MAF | `0.05` |
| Minimum sample size | `100` |
| COJO P-value threshold | `1e-5` |
| Locus window | `±100,000 bp` |
| SuSiE credible-set coverage | `0.95` |
| SuSiE maximum effects | `10` |

## Stage 1 outputs

Stage 1 writes all generated files under `data/finemap/female`:

```text
female_significant_snps.tsv
female_cojo_input.txt
bfile/merged_qc_numeric.bed
bfile/merged_qc_numeric.bim
bfile/merged_qc_numeric.fam
cojo/female_lifespan_cojo.jma.cojo
regions/chr<chromosome>_pos<position>_snps.tsv
```

The `.bed` and `.fam` files under `bfile` are symbolic links to the QC-filtered LD reference. The generated `.bim` contains numeric chromosome labels. Region files are replaced atomically, and regions no longer selected by COJO are removed.

## Stage 2 outputs

Stage 2 writes PLINK working files and SuSiE results under `data/finemap/female`:

```text
susie/chr<chromosome>_pos<position>.snplist
susie/chr<chromosome>_pos<position>.bed
susie/chr<chromosome>_pos<position>.bim
susie/chr<chromosome>_pos<position>.fam
susie/chr<chromosome>_pos<position>.ld
susie_results/chr<chromosome>_pos<position>_susie.tsv
```

Summary-statistic alleles are aligned to the PLINK BIM order before inference. Each result contains posterior inclusion probabilities and credible-set assignments. Obsolete SuSiE result files are removed after all current loci complete successfully.

## Tests

Run the pipeline test suite from the repository root:

```bash
python -m unittest tests/test_drosophila_female_lifespan_finemapping_pipeline.py -v
```
