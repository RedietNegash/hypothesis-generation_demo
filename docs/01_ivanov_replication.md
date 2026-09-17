# DGRP Female Lifespan: GWAS Replication and Heritability Partitioning

Replication of Ivanov et al. (2015) on DGRPool Study 18, extended with a
gene-level test (MAGMA) and cell-type heritability partitioning (LDSC-SEG).

Status: results below predate two pipeline corrections (section 7.1, 7.2) and
require regeneration before publication. See section 8.

Related: `02_ivanov_huang_merge.md` covers the cross-study phenotype merge,
which uses a different sample and phenotype scale.

## 1. Summary

| Analysis | Unit | Tests | Threshold | Significant |
|---|---|---|---|---|
| GWAS (plink2) | SNP | ~4.4M | 2.28e-8 | 2 at 5e-8, 0 at paper threshold |
| MAGMA | gene | 18,903 | 2.65e-6 | 0 |
| MAGMA gene-set | cell type | 163 | 3.07e-4 | 0 |
| LDSC-SEG | cell type | 163 | 3.07e-4 | 5 (10 at FDR 5%) |

SNP-level and gene-level analyses are null. Heritability partitioning is not:
five cell types carry significantly elevated per-SNP heritability. Section 5
reconciles these.

## 2. Data

| Field | Value |
|---|---|
| Phenotype | `S18_1537_F`, female lifespan in days, DGRPool Study 18 |
| Sample | 197 DGRP lines |
| Phenotype SD | 9.90 days |
| Genotypes | DGRP2 freeze 2, ~4.4M SNPs post-QC |
| Arms | 2L, 2R, 3L, 3R, 4, X |
| Cell-type annotations | 163 AFCA gene sets |

DGRP lines are inbred and effectively homozygous. N counts lines, not
individuals, and is not increasable by rearing more flies. This bounds all
downstream power (section 6).

## 3. Inputs from the paper

| Value | Source | Use |
|---|---|---|
| alpha = 2.28e-8 | `scripts/lifespan_power_analysis.py` | SNP significance threshold, adopted unchanged |
| R2 = 0.047 for all common variants | `scripts/gene_level_ftest_power.py` | power benchmark |

The threshold implies ~2.19M tests (0.05 / 2.28e-8). That count is inferred, not
sourced.

The paper's SNP-level findings are not recorded in this repository. Concordance
with the paper is therefore undetermined; results below stand on their own.

## 4. Results

### 4.1 GWAS

2 SNPs at p < 5e-8, 0 at 2.28e-8. Neither resolves to a credible set in
downstream fine-mapping, and neither is corroborated by 4.2.

Diagnostics: `results/lifespan_gwas_manhattan_qq.png` (Manhattan, QQ, lambda_GC),
`results/dgrp_pca_visualization.png` (PC1/PC2, scree).

### 4.2 MAGMA gene-level

18,903 genes, 1,356,952 of 1,965,595 SNPs mapped (69.0%). Zero significant under
Bonferroni, Holm-Bonferroni, or BH-FDR at 5%, 10%, 20%. Minimum p = 7.4e-5
against a 2.65e-6 threshold.

| Gene | NSNPS | ZSTAT | P |
|---|---|---|---|
| `INE-1{}6211` | 1 | 3.79 | 7.40e-5 |
| `INE-1{}5276` | 13 | 3.69 | 1.14e-4 |
| `Tdrd3` | 23 | 3.68 | 1.18e-4 |
| `Gcat` | 50 | 3.49 | 2.44e-4 |
| `snoRNA:Me18S-A28a` | 1 | 3.37 | 3.72e-4 |
| `Strica` | 17 | 3.33 | 4.29e-4 |

Two of the top six are single-SNP genes and two are INE-1 transposable element
annotations. The ranking is consistent with noise: real signal does not
concentrate in the smallest features of the annotation.

Gene-set test: 163 sets, 3,008 unique genes covered, zero significant.

### 4.3 LDSC-SEG

163 cell types, Bonferroni 0.05/163 = 3.07e-4. Five pass Bonferroni, ten pass
BH-FDR at 5%. All coefficients positive.

| Cell type | Coef | SE | P | q | Bonf |
|---|---|---|---|---|---|
| CNS surface-associated glial cell | 2.08e-5 | 5.18e-6 | 2.91e-5 | 0.0028 | yes |
| Female reproductive system | 2.44e-5 | 6.12e-6 | 3.47e-5 | 0.0028 | yes |
| Polar follicle cell | 1.87e-5 | 5.19e-6 | 1.59e-4 | 0.0087 | yes |
| Adult hindgut | 1.82e-5 | 5.26e-6 | 2.71e-4 | 0.0089 | yes |
| Enteroendocrine cell | 1.95e-5 | 5.65e-6 | 2.74e-4 | 0.0089 | yes |
| Pericerebral adult fat mass | 2.07e-5 | 6.20e-6 | 4.26e-4 | 0.0100 | no |
| Adult fat body, head | 2.25e-5 | 6.75e-6 | 4.29e-4 | 0.0100 | no |
| Epidermal cell, antimicrobial | 2.17e-5 | 6.94e-6 | 8.78e-4 | 0.0179 | no |
| Epithelial cell body | 1.82e-5 | 6.01e-6 | 1.20e-3 | 0.0206 | no |
| Oviduct | 1.91e-5 | 6.32e-6 | 1.26e-3 | 0.0206 | no |

Hits cluster into coherent tissue groups (glia and brain-adjacent fat,
reproductive, gut) rather than distributing arbitrarily across the 163.

Full table with corrections: `results/lifespan_female_cts_corrected.tsv`.

## 5. Interpretation

LDSC-SEG pools signal across every SNP near every gene in an annotation. MAGMA
requires one gene to clear 2.65e-6 alone. With common variants explaining 4.7%
of phenotypic variance in aggregate, the per-gene share is below detection at
N=197 while the per-annotation share is not. The null gene-level result and the
positive partitioning result are therefore consistent, not contradictory.

Supported: lifespan heritability in the DGRP is non-uniformly distributed and
concentrates in glial, reproductive, gut and fat annotations.

Not supported: any gene-level claim. LDSC-SEG resolves tissue, not genes, and
the analysis that would resolve genes returned nothing.

## 6. Power

Per-SNP, from `data/gwas/tmp/female_lifespan_power.tsv` (noncentral t,
alpha = 2.28e-8, SD = 9.90, N = 197):

| MAF | 5-day effect | 7.5-day | 10-day |
|---|---|---|---|
| 0.01 | 1.6e-7 | 5.8e-7 | 2.0e-6 |
| 0.05 | 3.1e-6 | 3.6e-5 | 3.1e-4 |
| 0.10 | 2.0e-5 | 4.0e-4 | - |

A 10-day effect is one phenotypic SD from a single locus. Power to detect it
remains below 0.1%.

Per-gene, from `scripts/gene_level_ftest_power.py`: MAGMA's gene test is an
F-test, so power is analytic. For a gene with k independent components,
ncp = N * R2 / (1 - R2) under noncentral F(k, N-k-1). The script computes the
R2 required for 80% power at 2.65e-6 per gene, using observed NPARAM values.
Compare against R2 = 0.047 for all common variants combined: any gene requiring
more than that individually is undetectable at this N.

The null gene-level result is a property of N=197, not of thresholding (the
Holm and FDR rows in 4.2 establish this) and not of implementation.

## 7. Methods

### 7.1 Genotype QC

```
plink2 --keep analysis_lines.keep --maf 0.01 --geno 0.05
```

`--keep` restricts QC to the 197 phenotyped lines. Frequency and missingness
must be computed on the analysis sample; otherwise a SNP can pass 1% MAF
panel-wide while being monomorphic among phenotyped lines, consuming a test
without contributing variance.

MAF 0.01 rather than the paper's 0.05: the LDSC annotation step (7.5) requires
SNP density, and at N=197 the frequency filter is not the binding constraint on
reliability (section 6 is).

### 7.2 Population structure

```
plink2 --indep-pairwise 200 50 0.2
plink2 --extract pca_prune.prune.in --pca 10
```

PC1 and PC2 enter the association model. LD pruning is required: the DGRP
segregates large cosmopolitan inversions, and unpruned components partly encode
inversion karyotype rather than ancestry, so conditioning on them removes real
signal from those regions. Pruning affects PC definition only; the GWAS tests
all QC'd SNPs.

Two components justified by scree flattening after PC2 and absence of discrete
clustering. Validated by lambda_GC.

### 7.3 Association

```
plink2 --linear hide-covar --covar-col-nums 3-4
```

Additive model per arm, `lifespan ~ SNP + PC1 + PC2`. Concatenated, Z = BETA/SE,
A2 derived from REF/ALT, munged via `munge_sumstats.py --signed-sumstats Z,0`.

### 7.4 MAGMA

Raw-genotype mode, `--annotate nonhuman`, 0 kb window. Multiple-testing
comparison across all result files: `scripts/magma_threshold_check.py`.

### 7.5 LDSC-SEG

`ldsc.py --h2-cts` over 163 cell types. LD scores at `--ld-wind-kb 1000` across
six arms. Two fly-specific substitutions:

- No baseline model exists for *Drosophila*. A genome-wide baseline was computed
  from all DGRP SNPs (`ldsc.py --l2`, no annotation) and passed as
  `--ref-ld-chr`.
- `--w-ld-chr` reuses that baseline. Human LDSC draws regression weights from a
  separate curated SNP list with no fly equivalent. Standard errors in 4.3 are
  consequently optimistic and the p-values are generous, not conservative.

## 8. Deviations from the paper

| Parameter | Paper | Here | Effect |
|---|---|---|---|
| MAF filter | 0.05 | 0.01 | ~4.4M vs ~2.19M tests |
| Ancestry covariates | study-specific | PC1 + PC2, LD-pruned | validated by lambda_GC |
| Wolbachia status | modelled | not modelled | see below |
| Inversion status | modelled | not modelled | see below |
| Gene-level test | not run | MAGMA | extension |
| Cell-type heritability | not run | LDSC-SEG | extension |

Wolbachia is a bacterial endosymbiont segregating across DGRP lines with a
direct lifespan effect. Inversions are large segregating haplotype blocks. Both
are standard DGRP covariates and neither is modelled here. LD pruning limits
inversion leakage into the PCs but does not substitute for modelling karyotype.
This is the largest outstanding gap.

## 9. Limitations

1. Results predate the QC restriction (7.1) and LD-pruned PCA (7.2). Both alter
   the SNP set and covariates. Section 4 requires regeneration; direction of
   change is not predictable since neither correction is a subset of prior
   behaviour.
2. Wolbachia and inversion covariates absent (section 8).
3. LDSC-SEG standard errors optimistic (7.5).
4. 0 kb MAGMA window excludes flanking regulatory variants. A 5 kb window was
   run separately; see `04_magma_framework.md`.
5. Three MAGMA runs exist (`magma_female_clean`, `clean2`, `clean3`) with
   differing gene counts. Only `clean3` is complete and cited. Delete or label
   the others.
6. Paper concordance undetermined (section 3).
7. `data/magma/`, `data/finemap/`, `data/gwas/` intermediates are absent from
   disk. Section 4.1 counts carry over from `ANALYSIS_PROGRESS.md`; 4.2 and 4.3
   were recomputed from `results/`.

## 10. Reproduction

```bash
jupyter nbconvert --execute notebooks/fly_ldsc_female_colab.ipynb

python scripts/run_magma_female_geneset.py              # 0 kb window
python scripts/run_magma_female_geneset.py --window 5   # 5 kb window
python scripts/magma_threshold_check.py

python scripts/lifespan_power_analysis.py
python scripts/gene_level_ftest_power.py
```

## 11. Artifacts

| Path | Contents |
|---|---|
| `results/lifespan_female_CellTypeSpecific.cell_type_results.txt` | raw h2-cts output, 163 rows |
| `results/lifespan_female_cts_corrected.tsv` | as above plus q-values and significance flags |
| `results/magma_female_clean3/dgrp_lifespan_female_gene.genes.out` | 18,903 gene results |
| `results/lifespan_female_cts_enrichment.png` | enrichment chart, Bonferroni line |
| `results/lifespan_gwas_manhattan_qq.png` | Manhattan, QQ, lambda_GC |
| `results/dgrp_pca_visualization.png` | PC1/PC2 scatter, scree |
| `data/gwas/tmp/female_lifespan_power.tsv` | per-SNP power grid |

Implementation: `notebooks/fly_ldsc_female_colab.ipynb`,
`scripts/run_magma_female_geneset.py`, `scripts/magma_threshold_check.py`,
`scripts/gene_level_ftest_power.py`, `scripts/lifespan_power_analysis.py`,
`scripts/lifespan_genotype_power.py`.
