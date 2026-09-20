# Fly Lifespan Analysis: What Was Tried, Methods, and Results

A chronological record of the analyses run on two fly lifespan datasets, the
methods used, what worked, and what did not. (For the condensed technical
reference, see `FINEMAPPING_SUMMARY.md`.)

---

## Datasets

| | DGRP female | Highfill DSPR |
|---|---|---|
| Source | DGRPool Study 18, trait `S18_1537_F` | Highfill et al. 2016, 805 pB DSPR RILs |
| Phenotype | Female lifespan (days) | Median lifespan (hours) |
| Sample size | **197 lines** | **803 RILs** |
| Panel type | Inbred wild-derived lines (DGRP2) | Recombinant inbred lines from 8 founders |
| Genotypes | ~4.4M SNPs (whole-genome sequenced) | 10,275 SNPs + 8-founder HMM probabilities |
| Original files | `data/gwas/lifespan_female_raw.tsv` | `data/RILSNP_R2.txt`, `data/Highfilletal(2016)805pBDSPRRILs.txt` |

Note the panels are fundamentally different: DGRP is a natural population
(individual SNP genotypes), DSPR is a designed cross (founder-haplotype
ancestry). This distinction turned out to matter a lot for which methods work.

---

## Methods tried, in order

### 1. GWAS (already run before this work)
- **DGRP female**: plink2 linear regression, `lifespan ~ SNP + PC1 + PC2`,
  4.4M SNPs. Result: **2 SNPs** at genome-wide significance (p<5e-8).
- **Highfill**: custom regression on a dosage matrix, with 8 PCs + Block as
  covariates, 9,125 SNPs. Result: **0 SNPs** at p<5e-8.
  - Note: the *uncorrected* version showed a strong hit at 3R:4292034
    (p=5.3e-9), but it vanished after PC correction -- a population-structure
    artifact, correctly removed.

### 2. MAGMA gene-based test
Aggregates all SNPs within each gene into a single gene-level p-value.

| | Genes tested | Bonferroni | BH-FDR<0.05 |
|---|---|---|---|
| DGRP female | 19,312 | **0** | **0** |
| Highfill | 2,225 | **2** (`mei-41`, `Fur2`) | **25** |

For DGRP female, also tested Holm-Bonferroni (FWER) and BH-FDR at 5/10/20%
per team suggestion -- **nothing passes at any threshold**. The closest gene
reaches q=0.21. This is a power limitation (N=197), not a thresholding choice.

### 3. MAGMA gene-set (cell-type) enrichment
163 cell-type gene sets tested in both datasets.
**Result: 0 significant sets in either dataset**, at any threshold.
Interpretation: the signal is a few strong single genes, not a broad
cell-type-specific program.

### 4. COJO (stepwise conditional selection of independent signals)
Identifies which associated SNPs are statistically independent of each other.

- **DGRP female**: 5 independent signals (used suggestive p<1e-5, since
  genome-wide 5e-8 leaves only 2 SNPs).
- **Highfill**: 3 independent signals -- X:16290404 (`mei-41`),
  3R:27570323, 2R:10101585.

### 5. SuSiE fine-mapping (credible sets)
Attempts to narrow each signal to a small set of likely causal variants.

- **DGRP female**: **1 credible set** (of 5 regions attempted) --
  2 variants, lead PIP 0.84. The other 4 regions produced no credible set.
- **Highfill**: produced 18 credible sets, but **these are not usable**
  (see "What did not work" below).

### 6. DSPRscan (founder-haplotype QTL mapping) -- Highfill only
The method actually designed for this experimental design: regresses
phenotype on the 8 founder-haplotype probabilities rather than SNP
genotypes. Required installing the real DSPR data package
(`DSPRqtlDataB`, ~3.2GB).

- 31 raw LOD peaks at threshold 6.8, which collapse to **4 distinct QTL
  regions** once peaks sharing a credible interval are merged
- 2 of the 4 are clean; 2 have degenerate founder estimates and were discarded
- Clean regions: **X:16.15-16.58Mb** (LOD 7.08) and **3R:4.04-5.02Mb**
  (LOD 9.35, strongest in the genome)

### 7. Variant-level resolution check within QTL intervals
The QTL intervals are wide (430kb and 980kb, containing 69 and 134
protein-coding genes), so the gene at the LOD peak is a guess, not an answer.
Checked whether finer resolution is achievable:
- SNPs inside the intervals are largely **independent** (median pairwise
  r2 = 0.03-0.04) -- so they *are* statistically separable
- **X region resolves**: three SNPs at p<2e-6 (`mei-41` 4.4e-7, `Fur2`
  5.1e-7, `mthl1` 1.7e-6), ~30x stronger than the next tier.
  **Narrows 69 genes to 3 candidates.**
- **3R region does not resolve**: no SNP passes even a region-level
  Bonferroni threshold (best p=1.7e-3), and the best SNP is in `CG45263`,
  not `hb`. The QTL is real but untagged by the sparse SNP panel
  (~1 SNP per 21kb). Needs denser genotyping, not a better algorithm.

---

## Results summary

### Confident findings
- **`mei-41`** (fly ortholog of human **ATR**, a DNA-damage checkpoint
  kinase) -- Bonferroni-significant in MAGMA, lead SNP of a COJO independent
  signal, and the top LOD peak of the founder-haplotype scan.
  **Three independent methods across two different genotype representations
  all converge on this gene.** Biologically coherent for an aging phenotype.
- **`Fur2`** -- Bonferroni-significant in MAGMA, adjacent to `mei-41` on X,
  and one of the 3 candidates the X interval narrows to.
- The X QTL narrows to **3 candidate genes**: `mei-41`, `Fur2`, `mthl1`
  (mutually correlated at r~0.8, so not separable from each other, but
  clearly distinguished from the other 66 genes in the interval).

### Weak / negative findings
- **`hb` (hunchback)** -- initially reported as a 3R candidate, but it has
  **no variant-level support**; it is only the gene at the LOD peak of a
  980kb interval containing 134 protein-coding genes. Treat as unsupported.
- **DGRP female produced no gene-level significance at all**, under any
  multiple-testing correction. Its single credible set maps to
  `lncRNA:CR44603` (a non-coding RNA, and the variants sit ~6.4kb *outside*
  it), which MAGMA independently ranks 2,942nd of 19,312 (p=0.16).
  The SNP-level and gene-level evidence disagree for this dataset.
- **No cell-type gene set is enriched** in either dataset.

---

## What did not work (and why)

**SuSiE on the Highfill dataset.** It ran and produced 18 credible sets,
but they are wrong: it assigned **PIP=0.0 to `mei-41`** -- the strongest and
most independently-validated signal in the dataset -- while giving PIP=1.0 to
several weaker neighbouring SNPs. SuSiE's own diagnostic also flagged
inconsistency between the summary statistics and the LD matrix.
Root cause: the LD reference had to be built from ~10% missing, discretized
allele-count-ratio calls rather than true sequencing genotypes -- too noisy
for SuSiE's joint Bayesian model. COJO's simpler stepwise approach tolerated
the same data; SuSiE did not. **This result should not be used.**

---

## Technical issues found and fixed along the way

These are worth recording because several would have silently produced
wrong answers:

1. **Missing dm3->dm6 liftover.** DSPR/Highfill positions are on the 2012
   dm3 assembly, but the gene annotation is dm6. Positions shift by up to
   **4.6Mb** -- without correcting this, QTL peaks map to entirely the wrong
   genomic region. Verified against a known-good reference file.
2. **Placeholder alleles in the existing MAGMA PLINK file.** 12 of 15 sampled
   SNPs had alleles that did not match the raw genotype calls (some were
   literally `0`). Fine for MAGMA (which only needs dosage) but unusable for
   COJO. Required rebuilding a genotype panel from raw read counts, with the
   dosage formula reverse-engineered and verified as
   `2 x countA1/(countA1+countA2)`.
3. **PLINK silently flipping alleles.** `--make-bed` reorders A1 to the minor
   allele by default, flipping **1,348 SNPs** away from the convention the
   GWAS effect sizes use. `--keep-allele-order` did *not* fix it; had to force
   assignment explicitly with `--a1-allele`.
4. **COJO's default collinearity threshold too lax for RIL data.** DSPR's
   long-range LD let 3 correlated SNPs (r=0.78-0.84) into one joint model,
   producing wildly inflated effect estimates (b~50 became bJ~300).
   Fixed with `--cojo-collinear 0.5`, which correctly collapses them to one.
5. **Monomorphic SNPs poisoning the LD matrix.** SNPs with no variation in a
   region have undefined correlation, which propagates NaN across every other
   SNP's row. Naive filtering wiped out entire regions; had to drop only the
   genuinely bad SNPs (NaN on the diagonal).
6. **Numerically degenerate variants and founder classes.** Near-monomorphic
   SNPs (and, in the founder analysis, founder groups with zero RILs) produce
   nonsensical estimates -- e.g. one founder effect estimated at -4,163,070
   hours with an SE of 17.5 million. Required explicit filtering and, for two
   DSPR regions, discarding the results.

---

## Scripts

| Script | Purpose |
|---|---|
| `notebooks/drosophila_female_lifespan_finemapping_pipeline.py` | DGRP female: COJO fine-mapping |
| `notebooks/fly_finemapping_susie_female.py` | DGRP female: SuSiE credible sets |
| `notebooks/fly_finemapping_gene_mapping_female.py` | DGRP female: credible set -> gene |
| `scripts/highfill_dspr_magma_gene.py` | Highfill: MAGMA gene + gene-set |
| `scripts/highfill_dspr_build_genotype_bfile.py` | Highfill: allele-consistent genotype panel |
| `scripts/highfill_dspr_finemapping.py` | Highfill: COJO fine-mapping |
| `scripts/highfill_dspr_susie.py` | Highfill: SuSiE (see caveat -- unreliable) |
| `scripts/highfill_dspr_qtl_scan.R` | Highfill: DSPR-native founder-haplotype QTL scan |

---

## Suggested next steps

1. **Denser genotyping for the 3R QTL** -- the strongest signal in the genome
   is currently untagged by any SNP. This is the clearest gap.
2. **Follow up `mei-41` biologically** -- it is the one triple-confirmed hit,
   and its human ortholog (ATR) has an established DNA-damage-repair role
   consistent with aging.
3. **Larger sample sizes for DGRP** -- N=197 is simply underpowered; no
   statistical correction rescues it.
4. **Feed the confident candidates into the hypothesis-generation pipeline**
   -- the original motivation for this fine-mapping work.
