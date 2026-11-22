# Run from the repository root: /usr/bin/Rscript rtd1_mtd/scripts/bulk_deseq2.R
project <- 'rtd1_mtd'
.libPaths(c(normalizePath(file.path(project, '.Rlib')), .libPaths()))
suppressPackageStartupMessages(library(DESeq2))
input <- file.path(project, 'data/processed')
sensitivity <- '--block-sensitivity' %in% commandArgs(trailingOnly=TRUE)
out <- file.path(project, if(sensitivity) 'results/bulk_sensitivity' else 'results/bulk')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
counts <- as.matrix(read.delim(file.path(input, 'human_counts_rounded.tsv'), row.names=1, check.names=FALSE))
samples <- read.delim(file.path(input, 'human_samples.tsv'), row.names=1, check.names=FALSE)
stopifnot(identical(colnames(counts), rownames(samples)), ncol(counts)==24,
          all(counts>=0), all(counts==round(counts)), !anyDuplicated(rownames(counts)))
samples$group <- factor(samples$group)
# GEO labels exp1–3 do not independently establish paired experimental blocks.
# Primary: use the eight treatment groups; do not infer pairing from label order.
# Sensitivity only: adjust for replicate labels as possible experimental blocks.
samples$replicate_label <- factor(samples$replicate_label)
design_formula <- if(sensitivity) ~replicate_label+group else ~group
keep <- rowSums(counts >= 10) >= 3
dds <- DESeqDataSetFromMatrix(counts[keep,], samples, design=design_formula)
dds <- DESeq(dds, parallel=FALSE)
genes <- read.delim(file.path(input, 'human_genes.tsv'))
summary <- list()
for (context in c('basal', 'LPS')) {
  for (dose in c('1','3','10')) {
    name <- paste0(context,'_',dose,'_vs_0')
    res <- as.data.frame(results(dds, contrast=c('group',paste0(context,'_',dose),paste0(context,'_0')), alpha=0.05))
    res$gene_id <- rownames(res)
    res <- merge(genes, res, by='gene_id', all.y=TRUE, sort=FALSE)
    res <- res[order(res$padj, na.last=TRUE),]
    write.table(res, file.path(out,paste0(name,'.tsv')), sep='\t', quote=FALSE, row.names=FALSE)
    summary[[name]] <- data.frame(contrast=name, tested=sum(!is.na(res$pvalue)),
                                 adjusted=sum(!is.na(res$padj)), significant=sum(res$padj<0.05,na.rm=TRUE),
                                 up=sum(res$padj<0.05 & res$log2FoldChange>0,na.rm=TRUE),
                                 down=sum(res$padj<0.05 & res$log2FoldChange<0,na.rm=TRUE))
  }
}
write.table(do.call(rbind,summary),file.path(out,'contrast_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pc <- prcomp(t(assay(vsd)))
pct <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
pdf(file.path(out,'sample_pca.pdf'),width=9,height=7)
palette <- setNames(hcl.colors(nlevels(samples$group),'Dark 3'),levels(samples$group))
plot(pc$x[,1:2],col=palette[samples$group],pch=19,xlab=paste0('PC1 (',pct[1],'%)'),ylab=paste0('PC2 (',pct[2],'%)'),main='Human THP-1: variance-stabilized counts')
text(pc$x[,1:2],labels=rownames(samples),pos=3,cex=0.6)
legend('topright',legend=levels(samples$group),col=palette,pch=19,cex=0.7)
dev.off()
write.table(data.frame(sample=rownames(samples),group=samples$group,PC1=pc$x[,1],PC2=pc$x[,2],size_factor=sizeFactors(dds)),file.path(out,'sample_qc.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(c(paste('Input genes:',nrow(counts)),paste('Retained genes:',sum(keep)),
 'Input: nearest-integer rounded Partek E/M estimated counts; not TMM-normalized expression.',
 paste('Filter: at least 10 counts in at least 3 samples. Design:',deparse(design_formula),'; eight groups, n=3 each.'),
 'Six within-LPS-stratum peptide versus vehicle contrasts. BH FDR 0.05 within each contrast.',
 if(sensitivity) 'Exploratory sensitivity only: replicate labels are candidate blocks, not verified matched biological replicates.' else 'Exploratory primary analysis; no MTD groups, no cross-species pooling, no inferred pairing.'),file.path(out,'methods.txt'))
print(do.call(rbind,summary))
