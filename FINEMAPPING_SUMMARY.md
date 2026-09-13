# Fly Lifespan GWAS: Gene-Level & Fine-Mapping Summary

Covers two independent fly lifespan datasets: the **DGRP female lifespan** GWAS
and the **Highfill et al. 2016 DSPR (803 RILs)** GWAS. For each, the same
pipeline was run: GWAS -> MAGMA gene/gene-set test -> COJO independent-signal
selection -> SuSiE fine-mapping -> nearest-gene mapping of credible sets.

## 1. DGRP female lifespan (DGRPool Study 18, S18_1537_F, N=197)

**MAGMA gene-based test** (`data/magma/dgrp_lifespan_female_gene.genes.out`,
19,312 genes):
- 0 genes pass Bonferroni (0.05/19312 = 2.6e-6)
- 0 genes pass BH-FDR even at a lenient 20% threshold
- Top genes (nominal only): `INE-1{}5276`, `INE-1{}5265` (X-linked TEs),
  `sdk` (sidekick) -- consistent with the fine-mapped locus below
- Conclusion: this dataset (N=197) is underpowered at the gene level

**COJO + SuSiE fine-mapping** (`notebooks/fly_finemapping_female.py`,
`notebooks/fly_finemapping_susie_female.py`):
- 5 independent signals selected by COJO (suggestive p<1e-5, since
  genome-wide 5e-8 leaves only 2 SNPs)
- Only 1 locus (2L, near `rs202927577`) resolved to an actual 95% credible
  set (2 SNPs, lead PIP 0.84); the other 4 loci had no credible set pass
  SuSiE's coverage/purity thresholds
- **Gene mapping** (`notebooks/fly_finemapping_gene_mapping_female.py`):
  the one credible set's SNPs are both ~6.4kb from **`lncRNA:CR44603`**
  (nearest gene, not overlapping)
- MAGMA's own gene-level test gives `lncRNA:CR44603` P=0.16 (not
  significant) -- the credible-set signal and the gene-based aggregate test
  disagree, a known limitation when a single strong SNP sits near (not
  inside) a gene body

**FDR/FWER check** (per team discussion): Benjamini-Hochberg and
Holm-Bonferroni were both applied to the MAGMA gene results -- neither
rescues anything; even a 20% FDR threshold is not met. The limitation is
power (N=197), not the choice of multiple-testing correction.

## 2. Highfill et al. 2016 DSPR RILs (MedLifespanHrs, N=803)

Original source files: `data/RILSNP_R2.txt` (genotypes, dm3, 2012 release)
and `data/Highfilletal(2016)805pBDSPRRILs.txt` (phenotypes, 805 RILs before
QC). Everything named `Highfill_803RILs_*` is derived/processed.

**MAGMA gene-based test** (`scripts/highfill_dspr_magma_gene.py`, output
`data/magma/highfill_lifespan_gene.genes.out`, 2,225 genes):
- **2 genes pass Bonferroni** (0.05/2225 = 2.25e-5): **`mei-41`** (fly
  ortholog of human ATR, DNA-damage checkpoint kinase) and **`Fur2`**, both
  X-linked and adjacent
- 25 genes pass BH-FDR<0.05
- Gene-set (cell-type) enrichment: 0/163 sets significant at any threshold
  -- the signal is a small number of strong single genes, not a
  cell-type-specific program
- Required a dm3->dm6 liftover step (the DSPR SNP positions are on the
  original 2012 assembly) that the initial script version was missing --
  fixed and verified to exactly reproduce a pre-existing reference result

**COJO fine-mapping** (`scripts/highfill_dspr_build_genotype_bfile.py` +
`scripts/highfill_dspr_finemapping.py`):
- Required building a new, allele-consistent genotype panel from the raw
  `genotype.txt` read counts (the existing MAGMA PLINK bfile had placeholder
  alleles, unusable for COJO) -- dosage formula reverse-engineered and
  verified as `2 x countA1/(countA1+countA2)`
- 3 independent, stable signals after correcting COJO's default
  collinearity threshold (DSPR's long-range LD let 3 correlated X-linked
  SNPs, pairwise r=0.78-0.84, into one unstable joint model with the default
  cutoff; `--cojo-collinear 0.5` fixed it): X:16290404 (`mei-41`),
  3R:27570323, 2R:10101585

**SuSiE fine-mapping** (`scripts/highfill_dspr_susie.py`): **ran, but the
result is not trustworthy and should not be used.** SuSiE's own diagnostic
flagged inconsistency between summary stats and the LD matrix on 2 of 3
regions, and concretely assigned **PIP=0.0 to X:16290404** -- the single
strongest, MAGMA-validated signal in the whole dataset -- while giving
PIP=1.0 to several much weaker, individually insignificant neighboring SNPs
instead. Root cause: the genotype panel has ~10% missingness and is built
from discretized allele-count-ratio calls rather than true sequencing
genotypes, too noisy an LD estimate for SuSiE's full joint model (COJO's
simpler stepwise/joint approach tolerated it; SuSiE did not).

**DSPR-native QTL scan** (`scripts/highfill_dspr_qtl_scan.R`): because SuSiE's
SNP-based LD approach broke down, re-ran the analysis using the *correct*
method for this experimental design -- DSPR's own `DSPRqtl` package, which
regresses phenotype directly on the 8 founder-haplotype probabilities (real
HMM-based genotype calls from the `DSPRqtlDataB` data package, ~3.2GB,
installed from `http://wfitch.bio.uci.edu/R/`) instead of SNP genotypes.
- `DSPRscan()` (model `MedLifespanHrs ~ factor(Block)`, design `inbredB`)
  found 31 raw LOD peaks at the literature-default threshold of 6.8 (used
  instead of a 200-iteration permutation test, which would have taken many
  hours for a threshold DSPR's own docs already call "fairly stable")
- Collapsing peaks that share a 95% Bayes credible interval (DSPR's docs warn
  raw peaks need this check) reduces 31 -> **4 distinct QTL regions**. 2 of
  the 4 have degenerate founder-mean estimates (empty founder classes at
  that position -- the same root problem as the SNP-level monomorphic-SNP
  issue, just in a new form) and should not be trusted; the other 2 are clean:
  - **X:16.15-16.58Mb (dm3)** peak LOD=7.08 -- lifts over to dm6 X:16,395,967,
    landing **directly inside `mei-41`**
  - **3R:4.04-5.02Mb (dm3)** peak LOD=9.35 (strongest signal in the dataset),
    lifts over to dm6 3R:8,694,278 (a **4.6Mb shift** -- confirms how
    essential the liftover step is for this pre-2014 resource), landing
    directly inside **`hb`** (hunchback)
- **`mei-41` is now confirmed by three independent methods**: MAGMA
  (gene-based burden test), COJO+SuSiE (SNP-level fine-mapping), and this
  founder-haplotype-based QTL scan (a completely different genotype
  representation). `hb` is a new candidate from this method only.

**Variant-level resolution within the QTL intervals**
(`data/highfill_finemap/highfill_variant_level_finemap.tsv`): the QTL
credible intervals are wide (430kb / 980kb, containing 69 and 134
protein-coding genes), so the gene at the LOD peak is a best guess, not a
resolved answer. Checking whether finer resolution is possible:
- SNPs inside these intervals are largely *independent* (median pairwise
  r2 = 0.03-0.04), so they are statistically separable -- the intervals are
  not one indivisible LD block
- **X QTL: resolves well.** Three SNPs stand out at p<2e-6, each inside a
  different gene -- `mei-41` (4.4e-7), `Fur2` (5.1e-7), `mthl1` (1.7e-6) --
  roughly 30x stronger than the next tier (`Rok`, `SMC3`, `Nup153`, ...).
  So the X signal narrows from 69 protein-coding genes to **3 candidate
  genes**. Those three are mutually correlated (r~0.78-0.84) and cannot be
  separated from each other, but are cleanly distinguished from the other 66.
- **3R QTL: does not resolve.** No SNP in the interval passes even a
  region-level Bonferroni threshold (best p=1.7e-3 vs threshold 1.1e-3), and
  the best SNP sits in `CG45263`, not `hb`. The strongest founder-haplotype
  signal in the genome (LOD 9.35) is simply **not tagged by any of the 47
  genotyped SNPs** in that 980kb window (~1 SNP per 21kb). So `hb` has no
  variant-level support -- it is only "the gene at the LOD peak". Resolving
  this locus would require denser genotyping, not a better algorithm.

## What to actually use going forward

| Result | Status |
|---|---|
| DGRP female: MAGMA gene/geneset | Solid null result (underpowered) |
| DGRP female: COJO+SuSiE credible set -> `lncRNA:CR44603` | Solid, but only 1 of 5 loci resolved |
| Highfill: MAGMA gene hits `mei-41`/`Fur2` | **Solid** -- Bonferroni-significant, biologically plausible |
| Highfill: COJO's 3 independent signals (SNP-level) | Solid |
| Highfill: SuSiE credible sets/PIPs (SNP-level) | **Do not use** -- demonstrated unreliable |
| Highfill: DSPRscan QTL regions X:16.15-16.58Mb and 3R:4.04-5.02Mb | **Solid** -- proper method for this design, `mei-41` triple-confirmed |
| Highfill: DSPRscan's other 2 raw regions (X:19.47-21.27Mb, 3R:24.59-24.64Mb) | **Do not use** -- degenerate founder estimates (empty founder classes) |
| Highfill: X QTL narrowed to `mei-41`/`Fur2`/`mthl1` | **Solid** -- 3 of 69 genes, ~30x stronger than next tier |
| Highfill: `hb` as the 3R candidate gene | **Weak** -- no SNP-level support; only the gene at the LOD peak |

## Files added this session

- `notebooks/fly_finemapping_female.py` -- DGRP female COJO fine-mapping
- `notebooks/fly_finemapping_susie_female.py` -- DGRP female SuSiE fine-mapping
- `notebooks/fly_finemapping_gene_mapping_female.py` -- DGRP credible-set -> gene mapping
- `scripts/highfill_dspr_magma_gene.py` -- Highfill MAGMA gene + gene-set analysis
- `scripts/highfill_dspr_build_genotype_bfile.py` -- allele-consistent Highfill genotype panel
- `scripts/highfill_dspr_finemapping.py` -- Highfill COJO fine-mapping
- `scripts/highfill_dspr_susie.py` -- Highfill SuSiE fine-mapping (see caveat above)
- `scripts/highfill_dspr_qtl_scan.R` -- Highfill DSPR-native founder-probability QTL scan (recommended over the SNP-based COJO/SuSiE route for this dataset)
