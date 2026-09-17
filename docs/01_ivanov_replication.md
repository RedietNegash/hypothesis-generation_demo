# 01: Replication of Ivanov et al. (2015), DGRP female lifespan GWAS

## What this document is

Three analyses were run on one dataset: a replication of the paper's GWAS, plus
two analyses the paper does not contain (a gene-level test and a cell-type
heritability test). This document explains what each one asks, what it returned,
and what the three together do and do not support.

**The short version.** The genome-wide scan and the gene-level test both came
back empty. The cell-type heritability test did not: five fly cell types carry
significantly more lifespan heritability than the genome average. That
combination is not a contradiction, and section 5 explains why.

A separate analysis merges this dataset with a second study's phenotypes. That
is a different method on a different sample and is documented in
`02_ivanov_huang_merge.md`. Nothing here uses it.

---

## 1. The starting point: what the paper established

Ivanov et al. ran a lifespan GWAS on the DGRP, a panel of inbred *Drosophila*
lines. Two numbers from that paper are used directly in this work:

| Value | Used here as |
|---|---|
| Significance threshold 2.28e-8 | the bar a SNP must clear, adopted unchanged |
| All common variants together explain about 4.7% of lifespan variance | the benchmark the power calculation is measured against |

Both are taken from the project's existing code
(`scripts/lifespan_power_analysis.py`, `scripts/gene_level_ftest_power.py`).
The threshold implies roughly 2.19 million SNPs tested (0.05 / 2.28e-8), but
that SNP count is a back-calculation, not something read off the paper.

**Open item: the paper's own SNP-level result is not recorded anywhere in this
repository.** Whether the replication below agrees or disagrees with the paper
cannot be stated until someone checks that. The results are reported here on
their own terms.

That 4.7% figure is the single most important number for interpreting
everything below. It says lifespan in the DGRP is highly polygenic: the genetic
signal is spread thinly across many variants, and no individual variant carries
much of it. A study designed to find individually large effects is therefore
looking for something the trait does not have.

---

## 2. The data, and the constraint that governs everything

| | |
|---|---|
| Phenotype | Female lifespan, DGRPool Study 18, trait `S18_1537_F` |
| Units | Days, untransformed |
| Sample | **197 DGRP lines** |
| Phenotype spread | mean lifespan varies with a standard deviation of 9.9 days |
| Genotypes | DGRP2 freeze 2, about 4.4 million SNPs after QC |
| Chromosome arms | 2L, 2R, 3L, 3R, 4, X |

The DGRP lines are inbred and effectively homozygous, so one line gives one
genotype and one phenotype value. **N = 197 counts lines, not flies.** Rearing
thousands of individual flies does not increase it. This is the binding
constraint on the whole analysis, and section 6 quantifies exactly how binding.

---

## 3. The three analyses, in plain terms

The three tests differ in what they treat as the unit of evidence: one SNP, one
gene, or one cell type's worth of genome. Reading left to right is reading from
most specific to most aggregated.

| | Unit tested | Question it answers | Result |
|---|---|---|---|
| GWAS | one SNP | Does any single variant shift lifespan on its own? | 2 SNPs at the conventional bar, none convincing |
| MAGMA | one gene | Do the SNPs in any gene, taken together, track lifespan? | nothing significant, out of 18,903 genes |
| LDSC-SEG | one cell type | Is heritability concentrated in the genes a given cell type uses? | **5 of 163 cell types significant** |

### 3.1 GWAS: testing one SNP at a time

For every SNP, fit a straight line predicting lifespan from genotype, while
holding population structure constant:

> lifespan = baseline + (effect of this SNP) + (correction for ancestry) + noise

Each SNP gets a p-value. With millions of SNPs tested, the bar is set at
2.28e-8 rather than the usual 0.05, because at 0.05 a scan of this size returns
tens of thousands of false hits by construction.

**Two corrections are built into this step.** Both are described in section 7.
Ancestry is handled by including principal components as covariates: without
them, a variant common in one subgroup of lines can look associated with
lifespan purely because that subgroup happens to be longer-lived.

**Result.** 2 SNPs reach p < 5e-8. Neither survives as a credible causal
variant in the follow-up fine-mapping (documented separately), and the
gene-level test below finds nothing corroborating them.

### 3.2 MAGMA: testing one gene at a time

A single SNP rarely carries much signal for a polygenic trait. MAGMA instead
asks whether the SNPs inside a gene are *collectively* more associated with
lifespan than chance would produce, yielding one p-value per gene. Testing
18,903 genes is a far smaller multiple-testing burden than testing 4.4 million
SNPs, so the bar drops to 0.05 / 18,903 = 2.65e-6.

**Result: nothing is significant.** Not under Bonferroni, not under
Holm-Bonferroni, and not under a false discovery rate of 5%, 10% or even 20%.
The strongest gene reaches p = 7.4e-5, well short of 2.65e-6.

The top of the ranking is worth reading, because it shows what a null result
looks like:

| Gene | SNPs in gene | P |
|---|---|---|
| `INE-1{}6211` | 1 | 7.40e-5 |
| `INE-1{}5276` | 13 | 1.14e-4 |
| `Tdrd3` | 23 | 1.18e-4 |
| `Gcat` | 50 | 2.44e-4 |
| `snoRNA:Me18S-A28a` | 1 | 3.72e-4 |
| `Strica` | 17 | 4.29e-4 |

Two of the top six are single-SNP "genes", and two are INE-1 transposable
element annotations rather than protein-coding genes. A real biological signal
does not preferentially land on the smallest and least gene-like features in
the annotation. This ranking is what random noise sorted by p-value looks like.

**Gene sets.** The same machinery was run on 163 cell-type gene sets, asking
whether the genes a cell type uses are collectively more associated than other
genes. Nothing significant.

### 3.3 LDSC-SEG: testing one cell type at a time

This test does not try to name a gene. It asks a statistical question about the
whole genome at once: **if you take all the genes a given cell type relies on,
do the SNPs near them account for more than their fair share of lifespan's
heritability?**

It can detect signal the other two cannot, because it never needs any single
SNP or gene to be individually significant. Thousands of variants each too weak
to name, but concentrated in the same set of genes, produce a detectable result
here and nothing at all in sections 3.1 and 3.2.

Each cell type gets a coefficient (positive means that cell type's regions carry
more heritability than average) and a p-value. 163 cell types were tested, so a
p-value of 0.05 is not evidence: about 8 cell types clear it by chance. The bar
is Bonferroni, 0.05 / 163 = 3.07e-4, with a false discovery rate reported
alongside.

**Result: 5 cell types pass Bonferroni, 10 pass FDR at 5%.**

| Cell type | Coefficient | Std. error | P | q (FDR) | Passes Bonferroni |
|---|---|---|---|---|---|
| CNS surface-associated glial cell | 2.08e-5 | 5.18e-6 | 2.91e-5 | 0.0028 | yes |
| Female reproductive system | 2.44e-5 | 6.12e-6 | 3.47e-5 | 0.0028 | yes |
| Polar follicle cell | 1.87e-5 | 5.19e-6 | 1.59e-4 | 0.0087 | yes |
| Adult hindgut | 1.82e-5 | 5.26e-6 | 2.71e-4 | 0.0089 | yes |
| Enteroendocrine cell | 1.95e-5 | 5.65e-6 | 2.74e-4 | 0.0089 | yes |
| Pericerebral adult fat mass | 2.07e-5 | 6.20e-6 | 4.26e-4 | 0.0100 | no |
| Adult fat body, head | 2.25e-5 | 6.75e-6 | 4.29e-4 | 0.0100 | no |
| Epidermal cell, antimicrobial response | 2.17e-5 | 6.94e-6 | 8.78e-4 | 0.0179 | no |
| Epithelial cell body | 1.82e-5 | 6.01e-6 | 1.20e-3 | 0.0206 | no |
| Oviduct | 1.91e-5 | 6.32e-6 | 1.26e-3 | 0.0206 | no |

Every coefficient is positive: each of these annotations absorbs more
heritability per SNP than the genome-wide average. The groupings are
biologically coherent rather than scattered (glia and brain-adjacent fat,
reproductive tissue, gut), which is what a real signal tends to look like and
what a random one usually does not.

Full table with both corrections: `results/lifespan_female_cts_corrected.tsv`.

---

## 4. Two adaptations that make this a fly analysis, not a human one

LDSC was built for human data and assumes two reference resources that do not
exist for *Drosophila*. Both substitutions are worth knowing when reading
section 3.3.

1. **No baseline model exists for the fly.** In human LDSC, enrichment is
   measured against a prebuilt baseline that accounts for generic genomic
   features. Here a baseline was computed from scratch across all DGRP SNPs and
   used in its place.
2. **The regression weights reuse that same baseline.** In human LDSC the
   weights come from a separate curated SNP list; no fly equivalent exists.
   The consequence is that the standard errors in the table above are
   **somewhat optimistic**, so the p-values should be read as slightly
   generous rather than conservative.

---

## 5. Reading the three results together

The apparent contradiction is: heritability is measurably structured by cell
type, yet not one gene is significant. Both can be true at once, and the reason
is arithmetic rather than biological.

LDSC-SEG pools a weak signal across every SNP near every gene in an annotation,
hundreds of genes at a time. MAGMA has to push a single gene past 2.65e-6 on
its own. For a trait where all common variants together explain only 4.7% of
variance, the per-gene share is minute, while the pooled share across a whole
cell type's gene repertoire is not.

**What this supports:** lifespan heritability in the DGRP is not uniformly
spread across the genome. It concentrates in genes used by glial, reproductive,
gut and fat tissue.

**What this does not support:** any claim about a specific gene. The cell-type
result names tissues, not genes, and the gene-level analysis that would name
them returned nothing.

---

## 6. Why the null results are the expected outcome

This section exists to distinguish "we found nothing" from "this design could
never have found anything". It is the second.

**Per-SNP power.** Given 197 lines, a phenotype standard deviation of 9.9 days,
and the paper's 2.28e-8 threshold, the probability of detecting a SNP that truly
does shift lifespan is:

| SNP frequency | True effect of 5 days | True effect of 10 days |
|---|---|---|
| 1% | 0.000016% | 0.0002% |
| 5% | 0.0003% | 0.03% |
| 10% | 0.002% | about 0.04% at 7.5 days |

A 10-day effect is a full standard deviation of lifespan from one locus, an
enormous effect for a polygenic trait. The chance of catching even that is far
below 1%. Source: `data/gwas/tmp/female_lifespan_power.tsv`.

**Per-gene power.** MAGMA's gene test is an F-test, so power can be computed
exactly rather than simulated. `scripts/gene_level_ftest_power.py` computes,
for each real gene, how much of lifespan variance that gene would have to
explain for an 80% chance of detection at 2.65e-6. The comparison is against
the paper's 4.7% figure for *all common variants combined*. Any gene requiring
more than that on its own is undetectable by the paper's own accounting.

**Therefore:** the empty gene-level result is a property of N = 197. It is not
a thresholding choice, which is exactly what the Holm and FDR rows in section
3.2 demonstrate, and not a pipeline defect. Only more lines would change it.

---

## 7. How each analysis was run

Reference material for reproducing or auditing the numbers above.

### 7.1 Genotype QC

```
--keep analysis_lines.keep   # restrict to the 197 phenotyped lines
--maf 0.01                   # drop SNPs with minor allele frequency below 1%
--geno 0.05                  # drop SNPs missing in more than 5% of lines
```

`--keep` matters: frequency and missingness must be computed on the same sample
the GWAS runs on. Without it, a SNP can clear the 1% frequency bar across the
full panel while being completely invariant among the phenotyped lines, and an
invariant predictor cannot explain anything but still consumes a test.

The paper used a 5% frequency filter. 1% is used here because the cell-type
analysis in section 3.3 needs SNP density more than it needs per-SNP
reliability, and because at N = 197 the frequency threshold is not what protects
the test (section 6 is).

### 7.2 Population structure

Merge all six arms, LD-prune (`--indep-pairwise 200 50 0.2`), compute 10
principal components on the pruned SNPs, use **PC1 and PC2** as covariates.

Pruning is not optional. The DGRP segregates several large cosmopolitan
inversions, which are long stretches of genome inherited as a block. Computed
without pruning, the leading components partly describe which inversions a line
carries rather than its overall ancestry, and correcting for those components
then removes real signal from inside those regions. Pruning affects only which
SNPs *define* the components; the GWAS still tests every QC'd SNP.

Two components are used because the scree plot flattens after PC2 and the
lines do not form discrete clusters. The check on that judgement is the genomic
inflation factor, printed by the notebook: a value near 1 means the correction
was sufficient.

### 7.3 Association testing

```
plink2 --linear hide-covar --covar-col-nums 3-4
```

Additive model, run per chromosome arm, results concatenated. Z = BETA / SE,
then passed through LDSC's `munge_sumstats.py`. `hide-covar` only suppresses
covariate rows from the output file; the covariates remain in the model.

### 7.4 Gene and gene-set analysis

MAGMA in raw-genotype mode, `--annotate nonhuman`, **0 kb window** (only SNPs
inside gene boundaries count). Multiple-testing comparisons across every result
file: `scripts/magma_threshold_check.py`.

### 7.5 Cell-type heritability

`ldsc.py --h2-cts` over 163 AFCA cell types, LD scores computed with
`--ld-wind-kb 1000` across all six arms, with the two fly adaptations from
section 4.

---

## 8. Differences from the paper

| | Paper | Here | Consequence |
|---|---|---|---|
| Frequency filter | 5% | 1% | more SNPs tested, about 4.4M vs 2.19M |
| Ancestry covariates | study-specific | PC1 + PC2, LD-pruned | checked via genomic inflation |
| Wolbachia infection status | accounted for | **not included** | see below |
| Chromosomal inversions | accounted for | **not included** | see below |
| Gene-level test | not run | MAGMA | extension |
| Cell-type heritability | not run | LDSC-SEG | extension |

**The Wolbachia and inversion gap is the most substantive difference.**
*Wolbachia* is a bacterial endosymbiont that infects some DGRP lines and affects
lifespan directly. Inversions are large genome blocks segregating in the panel.
Both are standard covariates in DGRP analyses, both affect the phenotype, and
neither is in this model. LD-pruning the PCA reduces inversion leakage into the
covariates but does not substitute for modelling inversion status. This is the
clearest item for a next iteration.

---

## 9. Limitations

1. **The results above predate two pipeline fixes.** QC is now restricted to
   phenotyped lines (7.1) and the PCA is now LD-pruned (7.2). Both change the
   SNP set and the covariates, so section 3's numbers need regenerating before
   being quoted as current. The direction of change is not predictable: neither
   fix is a strict subset of the old behaviour.
2. **Wolbachia and inversion covariates are missing** (section 8).
3. **LDSC-SEG standard errors are optimistic** because weights reuse the
   reference LD scores (section 4).
4. **The 0 kb MAGMA window** counts only SNPs inside gene bodies, missing
   regulatory variants nearby. A 5 kb window was run separately and is
   discussed in `04_magma_framework.md`.
5. **Three MAGMA runs exist** (`magma_female_clean`, `clean2`, `clean3`) with
   differing gene counts. `clean3` is the complete one and the only one cited
   here. The other two should be deleted or labelled.
6. **The paper's own SNP-level findings are not recorded in this repository**
   (section 1), so agreement or disagreement with it cannot yet be stated.
7. **Intermediate files under `data/magma/`, `data/finemap/` and `data/gwas/`
   are no longer on disk.** Section 3.1's SNP counts are carried over from
   `ANALYSIS_PROGRESS.md`. Sections 3.2 and 3.3 were recomputed from the
   surviving files in `results/`.

---

## 10. Reproducing

```bash
# Full pipeline, sections numbered as in the notebook
jupyter nbconvert --execute notebooks/fly_ldsc_female_colab.ipynb

# Gene and gene-set analysis
python scripts/run_magma_female_geneset.py              # 0 kb window
python scripts/run_magma_female_geneset.py --window 5   # 5 kb window

# Multiple-testing comparison across all result files
python scripts/magma_threshold_check.py

# Power
python scripts/lifespan_power_analysis.py
python scripts/gene_level_ftest_power.py
```

## 11. Files

| File | Contents |
|---|---|
| `results/lifespan_female_CellTypeSpecific.cell_type_results.txt` | raw cell-type output, 163 rows |
| `results/lifespan_female_cts_corrected.tsv` | the same, plus FDR q-values and both significance flags |
| `results/magma_female_clean3/dgrp_lifespan_female_gene.genes.out` | 18,903 gene-level results |
| `results/lifespan_female_cts_enrichment.png` | enrichment chart with the Bonferroni line |
| `results/lifespan_gwas_manhattan_qq.png` | Manhattan and QQ plots, with genomic inflation |
| `results/dgrp_pca_visualization.png` | PC1 vs PC2 scatter and scree plot |
| `data/gwas/tmp/female_lifespan_power.tsv` | per-SNP power grid |

**Implementation:** `notebooks/fly_ldsc_female_colab.ipynb` (end to end),
`scripts/run_magma_female_geneset.py`, `scripts/magma_threshold_check.py`,
`scripts/gene_level_ftest_power.py`, `scripts/lifespan_power_analysis.py`,
`scripts/lifespan_genotype_power.py`.
