# Female lifespan fine-mapping pipeline

The pipeline currently has three stages:

1. `fly_finemapping_female.py` selects independent GWAS signals with GCTA-COJO and extracts their loci.
2. `fly_finemapping_susie_female.py` fine-maps each locus with SuSiE-RSS.
3. `fly_finemapping_gene_mapping_female.py` maps selected variants to genes.

Only stage 1 is covered by the instructions below.

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

The Python environment must provide NumPy and pandas. GCTA must either be available as `gcta64` on `PATH` or supplied with `--gcta-bin`.

Record the Python and GCTA versions used for a run:

```bash
python --version
gcta64 --version
```

## Run stage 1

From the repository root, run:

```bash
python notebooks/fly_finemapping_female.py
```

For a GCTA executable outside `PATH`, run:

```bash
python notebooks/fly_finemapping_female.py \
  --gcta-bin /absolute/path/to/gcta64
```

An existing COJO result is reused by default. Use `--force` after changing inputs or analysis settings:

```bash
python notebooks/fly_finemapping_female.py \
  --gcta-bin /absolute/path/to/gcta64 \
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

## Tests

Run the stage 1 test suite from the repository root:

```bash
python -m unittest tests/test_fly_finemapping_female.py -v
```
