library(DSPRqtl)
library(DSPRqtlDataB)

base_dir <- "/mnt/hdd_1/rediet/fly-ldsc"

pheno <- read.table(
  file.path(base_dir, "data/Highfilletal(2016)805pBDSPRRILs.txt"),
  header = TRUE
)
names(pheno)[names(pheno) == "RIL"] <- "patRIL"

cat(sprintf("Loaded %d RILs\n", nrow(pheno)))

out_dir <- file.path(base_dir, "data/highfill_finemap/dspr_qtl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scan_file <- file.path(out_dir, "highfill_dsprscan.RData")
if (file.exists(scan_file)) {
  cat("Genome scan already exists, loading...\n")
  load(scan_file)
} else {
  cat("Running DSPRscan (genome scan against 8 founder probabilities)...\n")
  scan.results <- DSPRscan(
    MedLifespanHrs ~ factor(Block),
    design = "inbredB",
    phenotype.dat = pheno,
    id.col = "patRIL"
  )
  save(scan.results, file = scan_file)
}

cat("\nTop LOD scores:\n")
lod <- scan.results$LODscores
print(head(lod[order(-lod$LOD), ], 10))

# A full 200-iteration permutation test re-runs the whole genome scan 200
# times, which would take many hours. DSPR's own documentation says the
# literature default (6.8 for inbred designs) is "fairly stable" across
# phenotypes and recommended for initial data exploration, so use that
# directly instead. A proper permutation test can be run later (in parallel
# across the 28 available cores) if a data-specific threshold is needed.
threshold <- 6.8
cat(sprintf("\nUsing DSPR's literature-default threshold: %.1f\n", threshold))

cat("\nIdentifying QTL peaks...\n")
peaks <- DSPRpeaks(scan.results, threshold = threshold, LODdrop = 2)

cat(sprintf("\n%d significant QTL peak(s) found at threshold %.3f:\n", length(peaks), threshold))
for (i in seq_along(peaks)) {
  p <- peaks[[i]]
  bci <- p$CI$BCI  # p$CI is a list of two interval types (LODdrop, BCI);
                    # BCI is the 95% Bayes credible interval we actually want
  cat(sprintf(
    "\nPeak %d: chr %s, pos %d bp (LOD=%.2f), 95%% BCI: %s:%d-%d, %%var=%.2f, entropy=%.4f\n",
    i, p$peak$chr, p$peak$Ppos, p$peak$LOD,
    bci[1, "chr"], bci[1, "Ppos"], bci[2, "Ppos"],
    p$perct.var, p$entropy
  ))
  print(p$geno.means)
}

peaks_file <- file.path(out_dir, "highfill_dsprpeaks.RData")
save(peaks, file = peaks_file)
cat(sprintf("\nSaved: %s\n", peaks_file))
