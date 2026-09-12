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

## What to actually use going forward

| Result | Status |
|---|---|
| DGRP female: MAGMA gene/geneset | Solid null result (underpowered) |
| DGRP female: COJO+SuSiE credible set -> `lncRNA:CR44603` | Solid, but only 1 of 5 loci resolved |
| Highfill: MAGMA gene hits `mei-41`/`Fur2` | **Solid** -- Bonferroni-significant, biologically plausible |
| Highfill: COJO's 3 independent signals | Solid |
| Highfill: SuSiE credible sets/PIPs | **Do not use** -- demonstrated unreliable |

## Files added this session

- `notebooks/fly_finemapping_female.py` -- DGRP female COJO fine-mapping
- `notebooks/fly_finemapping_susie_female.py` -- DGRP female SuSiE fine-mapping
- `notebooks/fly_finemapping_gene_mapping_female.py` -- DGRP credible-set -> gene mapping
- `scripts/highfill_dspr_magma_gene.py` -- Highfill MAGMA gene + gene-set analysis
- `scripts/highfill_dspr_build_genotype_bfile.py` -- allele-consistent Highfill genotype panel
- `scripts/highfill_dspr_finemapping.py` -- Highfill COJO fine-mapping
- `scripts/highfill_dspr_susie.py` -- Highfill SuSiE fine-mapping (see caveat above)
