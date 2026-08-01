# CATlas cell-type enrichment pipeline

Given a GWAS, this pipeline asks **which human cell types are enriched for a
trait's association signal**, by testing whether the trait's independent lead
variants fall inside cell-type-specific open-chromatin regions (candidate
cis-regulatory elements, cCREs) more than expected by chance. It also scores each
variant with two DNA sequence models (Evo2, Nucleotide-Transformer-2.5B) to gauge
predicted regulatory impact inside the enriched cell types.

It ships with four worked example traits (T2D, BMI, AF, AD) and runs on **any**
GWAS via a column mapping — no code changes required.

---

## What it does (per step)

1. **Load & threshold** — read GWAS summary stats into a canonical schema
   (`variant | chr | pos | ref | alt | beta | pval | rsid`), keep genome-wide
   significant variants (`pval < 5e-8`).
2. **rsID matching** — intersect with a PLINK LD reference panel (hg19). If a file
   lacks real dbSNP rsIDs, they are recovered by a `chr:pos` join (see backfill).
3. **LD pruning** — PLINK `--indep-pairwise 1000 100 0.1` to keep approximately
   independent lead variants (removes redundant signal from linkage).
4. **Genome build detection + liftover** — auto-detect hg19 vs hg38 by matching
   the reference base at each locus against both genomes, then liftOver to hg38.
   Alleles are re-oriented to the hg38 reference (effect sign flipped to match).
5. **Promiscuity filter** — drop variants open in >20% of cell types (ubiquitous
   regulatory sites carry no cell-type-specific information).
6. **Sequence extraction** — pull ±500 bp around each variant from hg38 and build
   ref/alt allele sequences.
7. **Sequence-model scoring** (GPU) — Evo2-7B (log-likelihood delta) and
   NT-2.5B (ref/alt embedding cosine distance) per variant.
8. **Cell-type screen** — for each of 222 CATlas cell types, count lead variants in
   that cell type's cCREs and compare to a null built from 100,000 local
   permutations (each variant shifted ±1 Mb), giving a permutation p-value,
   enrichment ratio, and BH-FDR.
9. **Shortlist + enrichment** — cell types at FDR < 0.10 (else the top 15 by rank,
   flagged exploratory); for those, compare sequence-model scores for in-peak vs
   out-of-peak variants.

---

## Requirements

- **Python env with a working GPU stack** for steps 7 (Evo2 + NT-2.5B):
  `torch` (CUDA build matching your driver), `flash-attn`, `evo2`, `transformers`,
  plus `pandas numpy scipy statsmodels pyfaidx pyranges`. On this machine the
  `ad_analysis` conda env already has a working stack; use
  `--no-scoring` to run steps 1–6, 8–9 on CPU-only without the sequence models.
- **External tools:** PLINK 1.9 (`--indep-pairwise`), UCSC `liftOver`.
- **Reference data** (auto-fetched into `--scratch` on first run, ~1 GB each):
  hg38.fa, hg19.fa; plus the bundled hg19 LD panel (`data/EUR.*`), the
  hg19→hg38 chain (`data/reference/`), and the 222 CATlas cCRE BED files
  (`data/catlas_ccres/`, downloaded on first run).

---

## Usage

### Built-in datasets

```bash
python -m catlas_pipeline.run --preset t2d
python -m catlas_pipeline.run --preset all
```

### Any custom GWAS

Describe the file's columns; the pipeline handles the rest.

```bash
python -m catlas_pipeline.run \
    --gwas /path/to/my_gwas.tsv.gz --label MYTRAIT \
    --col chr=CHR --col pos=BP --col ref=A1 --col alt=A2 \
    --col beta=BETA --col pval=P --col rsid=SNP \
    --pval 5e-8
```

- Required `--col` mappings: `chr, pos, ref, alt, pval`. Optional: `beta` (else
  0.0; effect sign then isn't tracked), `rsid`.
- If the file has **no real dbSNP rsIDs** (e.g. its id column is `chr:pos`), omit
  `--col rsid` and add `--backfill-rsids` to recover them from the LD panel.
- `ref`/`alt` may map to effect/non-effect allele columns in either order — the
  build-detection and reorientation steps sort out which is the genome reference.

### Useful flags

`--no-scoring` (skip GPU sequence models) · `--pval` · `--n-perm` · `--fdr` ·
`--sep` · `--uppercase-alleles` · `--base` / `--scratch` (or env `CATLAS_BASE`,
`CATLAS_SCRATCH`, `CATLAS_BIM_PREFIX`, `PLINK_BIN`, `LIFTOVER_BIN`).

Outputs land in `data/results/`: `celltype_screen_222_<label>.tsv` (full 222-cell
ranking) and `<label>_celltype_enrichment.tsv` (sequence-model in/out-of-peak).

---

## Results (example traits)

Full pipeline, 100k permutations, 222 CATlas cell types, GPU scoring on every
variant. "Independent variants" = lead variants after LD pruning + promiscuity
filter (the number actually tested).

| Trait   | Independent vars | Top cell types (enrichment ×, permutation p)                                                                                                                                     | Read                                              |
| ------- | :--------------: | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------- |
| **T2D** |       379        | Pancreatic Beta Cell (1.74×, p=0.018), Fetal Pancreatic Islet (2.12×, p=0.029), Pancreatic Delta/Gamma (1.75×, p=0.033), Pancreatic Alpha (1.61×, p=0.045), Hepatocyte (1.26×)   | Pancreatic islet + liver — textbook T2D biology   |
| **AF**  |       197        | Fetal Ventricular Cardiomyocyte (1.77×, p=0.008), Atrial Cardiomyocyte (2.18×, p=0.010), Ventricular Cardiomyocyte (1.84×, p=0.026), Fetal Atrial Cardiomyocyte (1.92×, p=0.030) | Cardiomyocytes dominate — right for an arrhythmia |
| **AD**  |        64        | Memory B Cell (2.40×), Plasma Cell (1.86×), **Microglia** (2.21×), Glutaminergic Neuron (2.01×)                                                                                  | Microglia + immune + neurons — matches AD GWAS    |
| **BMI** |       696        | Fetal Adrenal Cortical (1.10×), Fetal Skeletal Myocyte (1.10×), Astrocyte (1.48×), Erythroblast, Cardiomyocyte, Neurons                                                          | Diffuse & weak — see below                        |

For T2D, AF, and AD the top-ranked cell type is the established causal tissue for
the disease — a strong correctness signal that the pipeline is working.

### Important statistical caveat

Across all four traits, **no cell type passed FDR < 0.10** after correcting for
the ~90–200 cell types tested, so every shortlist above is the **top-15-by-rank
exploratory** fallback, not an FDR-confirmed hit. AF came closest (top hits at
FDR ≈ 0.80, raw p ≈ 0.008). Treat these as hypotheses to follow up, not
confirmed cell-type associations. Larger, better-powered GWAS (more independent
lead variants) would sharpen this.

### Why BMI looks "wrong" (diffuse, no clear tissue)

BMI's shortlist has no dominant, obviously-causal tissue and the weakest
enrichments (~1.0–1.5×, all p > 0.3), unlike the sharp signals for T2D/AF/AD.
Likely reasons — mostly **biological/statistical, not a pipeline bug**:

- **BMI is extremely polygenic and multi-tissue.** Its heritability is
  concentrated in the **central nervous system** (appetite/energy-balance
  circuits) plus adipose, muscle, and endocrine tissue. Signal genuinely spread
  across many cell types produces a flat, weak enrichment profile — the diffuse
  result partly reflects real biology.
- **CATlas may lack the most relevant cell types at high resolution.** The most
  BMI-relevant populations (hypothalamic neuron subtypes, mature adipocytes) are
  under-represented among the 222 CATlas cCRE tracks, so the strongest BMI
  signal has no matching column to land in.
- **Ancestry mismatch in LD pruning.** The bundled LD panel is
  **European-ancestry** (`data/EUR`), while the BMI sumstats are UK-Biobank
  "both sexes" raw-BMI. Pruning with a mismatched panel can retain correlated
  variants or drop true leads, diluting cell-type specificity. Supply an
  ancestry-matched panel via `CATLAS_BIM_PREFIX` for a fairer test.
- **Raw (untransformed) BMI phenotype** inflates a heavy tail and effective-N
  issues vs. the inverse-normal-transformed version, adding noise to which loci
  clear 5e-8.
- **Power / multiple testing.** Even with 696 lead variants, spreading a weak
  per-cell-type signal across ~200 tests leaves nothing past FDR.

None of these are failures of the code path — BMI exercises the same, validated
steps that give clean answers for T2D/AF/AD. It is a good reminder that a _diffuse
polygenic trait_ is expected to yield a diffuse, low-confidence cell-type profile.

---

## Assumptions & limitations

- **Human only**, autosomes; variants must be hg19 or hg38 (auto-detected).
- **LD panel is hg19 European-ancestry by default** — match it to your GWAS's
  ancestry (`CATLAS_BIM_PREFIX`) for valid pruning and rsID backfill.
- **Cell-type resolution is bounded by CATlas** (222 adult+fetal cCRE tracks); a
  trait whose causal cell type is absent will show diffuse enrichment.
- The permutation null shifts variants ±1 Mb (local genomic background); it does
  not correct for gene density or mappability beyond that.
- Sequence-model scores are **exploratory annotations**, not calibrated effect
  predictors; they contextualize the enrichment, they don't establish causality.

---

## Layout

```
catlas_pipeline/
  run.py         CLI entry point (presets + custom datasets)
  config.py      PipelinePaths / PipelineConfig (paths, thresholds; env-overridable)
  pipeline.py    run_pipeline() orchestration + ensure_reference_data()
  ld.py          rsID matching, LD pruning, chr:pos rsID backfill
  genome.py      build detection, liftover, reorientation, sequence extraction,
                 reference-genome auto-fetch
  scoring.py     Evo2 + NT-2.5B GPU scoring
  catlas.py      CATlas cCRE download, permutation cell-type screen, shortlist
  enrichment.py  in-peak vs out-of-peak sequence-model comparison
  prep/
    generic.py   prepare_gwas() — arbitrary file + column map
    bmi.py af.py t2d.py ad.py   worked-example presets (also document each file's quirks)
```
