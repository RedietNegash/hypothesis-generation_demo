#!/bin/bash
for CHR in chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do
    COUNT=$(ls /mnt/hdd_1/rediet/deltaSVM/runs/$CHR/tmp/*.merge.gkm.tsv 2>/dev/null | wc -l)
    PARQUET=$(ls /mnt/hdd_1/rediet/deltaSVM/runs/$CHR/out/*_deltasvm.parquet 2>/dev/null && echo "done" || echo "not done")
    echo "$CHR: $COUNT/94 TFs | $PARQUET"
done