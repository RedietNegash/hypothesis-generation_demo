# System R >=4.6; package library is confined to this project.
stopifnot(getRversion() >= '4.6')
dir.create('rtd1_mtd/.Rlib', recursive=TRUE, showWarnings=FALSE)
.libPaths(c(normalizePath('rtd1_mtd/.Rlib'), .libPaths()))
# Rebuild these compiled packages to avoid pre-existing libraries from older R.
install.packages(c('Rcpp','RcppArmadillo','rlang','cli','vctrs','glue','ggplot2'),
                 lib='rtd1_mtd/.Rlib', repos='https://cloud.r-project.org', Ncpus=2)
install.packages('DESeq2', lib='rtd1_mtd/.Rlib',
                 repos=c('https://bioconductor.org/packages/3.23/bioc','https://cloud.r-project.org'), Ncpus=2)
stopifnot(requireNamespace('DESeq2', quietly=TRUE))
