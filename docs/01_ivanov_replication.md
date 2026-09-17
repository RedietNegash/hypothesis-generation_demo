# 01 — Replication of Ivanov et al. (2015): DGRP female lifespan GWAS

**Scope.** A single-dataset replication of the longevity GWAS described in
*Longevity GWAS Using the Drosophila Genetic Reference Panel*, followed by two
extensions the paper does not contain: a MAGMA gene/gene-set analysis and an
LDSC-SEG cell-type heritability analysis on the same summary statistics.

The cross-study merge with the Huang female data is a **different method with a
different phenotype scale and sample size** and is documented separately in
`02_ivanov_huang_merge.md`. Nothing in this document uses the merged phenotype.

**Implementation.** `notebooks/fly_ldsc_female_colab.ipynb` (end-to-end),
`scripts/run_magma_female_geneset.py`, `scripts/magma_threshold_check.py`,
`scripts/gene_level_ftest_power.py`, `scripts/lifespan_power_analysis.py`,
`scripts/lifespan_genotype_power.py`.

---

## 1. Data

| | |
|---|---|
| Phenotype | Female lifespan, DGRPool Study 18, trait `S18_1537_F` |
| Units | Days (raw, untransformed) |
| Lines with phenotype ∩ genotype | **197** |
| Genotypes | DGRP2 freeze 2, one PLINK fileset per chromosome arm |
| Arms analysed | 2L, 2R, 3L, 3R, 4, X (autosomes + X; no heterochromatin) |
| Phenotype SD | 9.90 days |

The DGRP lines are inbred and effectively homozygous, so each line contributes
one genotype and one phenotype mean. `N` is a count of **lines**, not flies —
this is what makes the panel small and is the single dominant constraint on
everything below.

## 2. Statistical methods

### 2.1 Genotype QC
Per chromosome arm, with plink2:

```
--keep analysis_lines.keep   # restrict to the phenotyped lines
--maf 0.01                   # drop SNPs with minor allele frequency < 1%
--geno 0.05                  # drop SNPs missing in > 5% of lines
```

`--keep` matters: allele frequency and missingness must be computed on the
sample the GWAS actually runs on. Without it a SNP can clear 1% MAF across the
full DGRP2 panel while being monomorphic among the phenotyped lines, and a
monomorphic predictor contributes nothing but still consumes a test.

Rationale for MAF 1% rather than the more common 5%: at N=197 a 5% MAF
threshold still leaves only ~10 minor-allele lines, so the threshold is not
what protects the test — the effect-size requirement is (§5). Keeping 1%
preserves rare-variant coverage for the LDSC annotation step, which needs SNP
density more than it needs per-SNP reliability. **This is a deviation from the
paper**, which used a 5% cut-off on ~2.19M SNPs.

### 2.2 Population structure
All six QC'd arms are merged into one genome-wide fileset, LD-pruned
(`--indep-pairwise 200 50 0.2`), and the top 10 principal components are
computed on the pruned SNPs (`plink2 --pca 10 --extract`). **PC1 and PC2 only**
are carried into the association model; the GWAS itself still tests every QC'd
SNP, pruning only decides which SNPs define the components.

Pruning is not optional here. The DGRP segregates several large cosmopolitan
inversions (*In(2L)t*, *In(3R)Mo* and others) which are long, high-LD blocks;
computed on unpruned genotypes the leading PCs partly describe inversion
karyotype rather than genome-wide ancestry, and conditioning on them then
removes real signal from inside those regions.

Justification for stopping at two PCs: the scree plot
(`results/dgrp_pca_visualization.png`) shows variance explained flattening
after PC2, and the PC1/PC2 scatter shows the DGRP lines distributed uniformly
rather than in discrete clusters — this panel has no strong stratification to
remove. The empirical check is the genomic inflation factor (§3.1), which is
what would expose insufficient correction.

### 2.3 Association model
Per chromosome arm:

```
plink2 --linear hide-covar --covar-col-nums 3-4
```

which fits, for each SNP *j*:

> lifespan_i = β₀ + β_j · g_ij + γ₁·PC1_i + γ₂·PC2_i + ε_i

with `g_ij ∈ {0, 2}` in practice (inbred lines), additive coding, and ε assumed
i.i.d. normal. `hide-covar` suppresses covariate rows from the output; the
covariates are still in the model. Per-arm results are concatenated,
Z = BETA/SE is computed, A2 is joined from the `.bim`, and the table is passed
through LDSC's `munge_sumstats.py` with `--signed-sumstats Z,0`.

### 2.4 Significance thresholds
Three thresholds are reported side by side, deliberately:

| Threshold | Value | Source |
|---|---|---|
| Bonferroni, paper | 2.28e-8 | 0.05 / 2,193,000 SNPs tested by Ivanov et al. |
| Genome-wide, conventional | 5e-8 | field standard |
| Suggestive | 1e-5 | used downstream to seed COJO, which has nothing to condition on at 5e-8 |

The paper's own threshold is reused verbatim so the two studies are judged on
identical terms rather than on whichever cut-off happens to flatter this run.

### 2.5 Gene-level and gene-set analysis (extension, not in the paper)
MAGMA in raw-genotype mode: `--annotate nonhuman` with a **0 kb window**,
gene-based test on the QC'd PLINK fileset, then a competitive gene-set test
over 163 AFCA cell-type gene sets. Multiple testing handled in
`scripts/magma_threshold_check.py`: Bonferroni, Holm–Bonferroni (FWER), and
Benjamini–Hochberg at q = 0.05, 0.10 and 0.20.

### 2.6 Cell-type heritability (extension, not in the paper)
LDSC-SEG (`ldsc.py --h2-cts`) over the same 163 cell types. Two fly-specific
adaptations were unavoidable:

- **No baseline model exists for Drosophila.** A genome-wide baseline LD score
  file was built from all DGRP SNPs with `ldsc.py --l2` and no annotation, and
  passed as `--ref-ld-chr`.
- **The same baseline is used for `--w-ld-chr`.** In human LDSC the regression
  weights come from a separate HapMap3-restricted set; no equivalent curated
  SNP list exists for the fly, so weights and reference are the same file. This
  makes the standard errors mildly optimistic and is a known limitation.

LD scores were computed with `--ld-wind-kb 1000` across all six arms.

### 2.7 Correction across cell types
163 cell types are tested in one `h2-cts` run, so nominal p < 0.05 is not a
significance claim. Primary threshold is **Bonferroni, 0.05/163 = 3.07e-4**,
with **BH-FDR** reported alongside. Corrected table:
`results/lifespan_female_cts_corrected.tsv`.

## 3. Results

### 3.1 Single-SNP GWAS
- ~4.4M SNPs tested after QC; **2 SNPs** reach p < 5e-8.
- The gene-level analysis (§3.2) finds nothing, and no SNP-level signal
  survives into a credible set except the one described in the fine-mapping
  document — so the replication outcome is a **negative result that matches
  the paper's own negative result**, not a contradiction of it.
- λ_GC and the Manhattan/QQ plots are produced in the notebook
  (`results/lifespan_gwas_manhattan_qq.png`); λ near 1 is what confirms the
  two-PC correction is sufficient.

### 3.2 MAGMA gene and gene-set
Run: `results/magma_female_clean3/` (0 kb window).

| | |
|---|---|
| Genes tested | **18,903** |
| SNPs mapped to genes | 1,356,952 of 1,965,595 in the fileset (69.0%) |
| Bonferroni threshold | 0.05 / 18,903 = 2.65e-6 |
| Genes passing Bonferroni | **0** |
| Genes passing Holm–Bonferroni | **0** |
| Genes passing BH-FDR at 5 / 10 / 20% | **0 / 0 / 0** |
| Smallest gene p-value | 7.40e-5 (`INE-1{}6211`) |

Top genes by p-value (nominal only, none significant):

| Gene | Chr | NSNPS | ZSTAT | P |
|---|---|---|---|---|
| `INE-1{}6211` | 2 | 1 | 3.79 | 7.40e-5 |
| `INE-1{}5276` | X | 13 | 3.69 | 1.14e-4 |
| `Tdrd3` | 3 | 23 | 3.68 | 1.18e-4 |
| `Gcat` | 3 | 50 | 3.49 | 2.44e-4 |
| `snoRNA:Me18S-A28a` | 2 | 1 | 3.37 | 3.72e-4 |
| `Strica` | 2 | 17 | 3.33 | 4.29e-4 |

Two of the top six are single-SNP genes (`NSNPS = 1`) and two are INE-1
transposable-element annotations — both patterns are what a null distribution
produces, not evidence.

**Gene sets:** 163 sets read, 3,008 unique genes covered. **0 sets significant
at any threshold.** MAGMA truncates Z-scores 3 points below zero and 6 SDs
above the mean before fitting, as logged.

### 3.3 Cell-type heritability (LDSC-SEG)
163 cell types; **5 pass Bonferroni, 10 pass BH-FDR < 0.05.**

| Cell type | Coefficient | SE | P | q | Bonferroni |
|---|---|---|---|---|---|
| CNS surface-associated glial cell | 2.08e-5 | 5.18e-6 | 2.91e-5 | 0.0028 | ✓ |
| Female reproductive system | 2.44e-5 | 6.12e-6 | 3.47e-5 | 0.0028 | ✓ |
| Polar follicle cell | 1.87e-5 | 5.19e-6 | 1.59e-4 | 0.0087 | ✓ |
| Adult hindgut | 1.82e-5 | 5.26e-6 | 2.71e-4 | 0.0089 | ✓ |
| Enteroendocrine cell | 1.95e-5 | 5.65e-6 | 2.74e-4 | 0.0089 | ✓ |
| Pericerebral adult fat mass | 2.07e-5 | 6.20e-6 | 4.26e-4 | 0.0100 | — |
| Adult fat body (head) | 2.25e-5 | 6.75e-6 | 4.29e-4 | 0.0100 | — |
| Epidermal cell, antimicrobial response | 2.17e-5 | 6.94e-6 | 8.78e-4 | 0.0179 | — |
| Epithelial cell body | 1.82e-5 | 6.01e-6 | 1.20e-3 | 0.0206 | — |
| Oviduct | 1.91e-5 | 6.32e-6 | 1.26e-3 | 0.0206 | — |

All coefficients are positive, i.e. each of these annotations absorbs more
per-SNP heritability than the genome-wide average.

### 3.4 The central tension in this replication
**The heritability is detectably structured by cell type, but no individual
gene is significant.** That is internally consistent rather than contradictory:
LDSC-SEG aggregates a weak polygenic signal across an entire annotation, while
MAGMA must resolve a single gene above a 2.65e-6 bar. At N=197 the second is
out of reach even when the first is not. The quantitative version of that
argument is §5.

## 4. Deviations from the paper

| | Paper | Here | Why |
|---|---|---|---|
| MAF filter | 5% | 1% | preserves SNP density for LDSC annotations; does not change power at this N |
| SNPs tested | ~2.19M | ~4.4M | consequence of the MAF filter |
| Covariates | study-specific | PC1 + PC2 | scree/scatter show no further structure |
| Gene-level test | — | MAGMA | extension |
| Cell-type heritability | — | LDSC-SEG | extension |
| Wolbachia / inversion covariates | accounted for | **not included** | known gap, see §6 |

## 5. Power — why the negative results are expected

Two independent calculations, neither of which needs a simulation.

**Per-SNP power** (`scripts/lifespan_power_analysis.py`, noncentral *t*, α =
2.28e-8, SD = 9.90 days, N = 197). Output: `data/gwas/tmp/female_lifespan_power.tsv`.

| MAF | 5-day effect | 10-day effect |
|---|---|---|
| 0.01 | 1.6e-7 | 2.0e-6 |
| 0.05 | 3.1e-6 | 3.1e-4 |
| 0.10 | 2.0e-5 | ~4e-4 (at 7.5 days) |

A 10-day shift in lifespan is a **one-SD** effect from a single locus. Power to
detect even that is well under 1%.

**Per-gene power** (`scripts/gene_level_ftest_power.py`). MAGMA's gene test is
a PC-regression F-test, so power is analytic: with *k* = NPARAM independent
components and N samples, the statistic is F(k, N−k−1) under the null and
noncentral F with **ncp = N·R²/(1−R²)** under an effect of size R². The script
computes, per real gene, the R² needed for 80% power at α = 0.05/18,903.

The benchmark is the paper's own estimate that **all common variants together
explain ~4.7% of lifespan variance**. Any single gene needing an R² above that
to be detectable is, by the paper's own accounting, undetectable here.

**Conclusion: the null gene-level result is a property of N=197, not of the
pipeline or the choice of multiple-testing correction.** No threshold rescues
it — that is exactly what §3.2's Holm and BH rows demonstrate.

## 6. Known limitations

1. **No Wolbachia infection status or chromosomal inversion covariates.** Both
   are standard in DGRP analyses and both are absent here. This is the most
   substantive gap relative to the paper.
2. **`--w-ld-chr` reuses the reference LD scores** (§2.6); LDSC-SEG standard
   errors are mildly optimistic as a result.
3. **0 kb MAGMA window** misses regulatory variants outside gene bodies. A
   ±5 kb window was run separately and is discussed in `04_magma_framework.md`.
4. **Three MAGMA runs exist** (`magma_female_clean`, `clean2`, `clean3`) with
   differing gene counts. `clean3` is the complete one and the only one cited
   here; the earlier two should be deleted or labelled.
5. **The results in this document predate two pipeline fixes.** QC is now
   restricted to the phenotyped lines (§2.1) and PCA is now LD-pruned (§2.2);
   both change the SNP set and the covariates, so every number in §3 needs
   regenerating before it is quoted as current. The direction of change is not
   predictable in advance — neither fix is a strict subset of the old
   behaviour.
6. **Intermediate files under `data/magma/`, `data/finemap/` and
   `data/gwas/` are no longer on disk.** GWAS-level counts in §3.1 are carried
   over from `ANALYSIS_PROGRESS.md`; everything in §3.2 and §3.3 was recomputed
   from surviving files in `results/` while writing this document.

## 7. Reproducing

```bash
# End-to-end pipeline (sections numbered as in the notebook)
jupyter nbconvert --execute notebooks/fly_ldsc_female_colab.ipynb

# Gene / gene-set analysis
python scripts/run_magma_female_geneset.py              # 0 kb window
python scripts/run_magma_female_geneset.py --window 5   # +/-5 kb window

# Multiple-testing report across all result files
python scripts/magma_threshold_check.py

# Power
python scripts/lifespan_power_analysis.py
python scripts/gene_level_ftest_power.py
```

## 8. Output files

| File | Contents |
|---|---|
| `results/lifespan_female_CellTypeSpecific.cell_type_results.txt` | raw `h2-cts` output, 163 cell types |
| `results/lifespan_female_cts_corrected.tsv` | as above + BH q-values and both significance flags |
| `results/magma_female_clean3/dgrp_lifespan_female_gene.genes.out` | 18,903 gene-level results |
| `results/lifespan_female_cts_enrichment.png` | enrichment bar chart with the Bonferroni line |
| `results/lifespan_gwas_manhattan_qq.png` | Manhattan + QQ with λ_GC |
| `results/dgrp_pca_visualization.png` | PC1/PC2 scatter + scree |
| `data/gwas/tmp/female_lifespan_power.tsv` | per-SNP power grid |
