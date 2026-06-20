#!/bin/bash

CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
BASE_DIR="/mnt/hdd_1/rediet/deltaSVM"

for CHROM in $CHROMS; do
    WORK_DIR="$BASE_DIR/runs/$CHROM"
    mkdir -p "$WORK_DIR"

    cp -r "$BASE_DIR/scripts" "$WORK_DIR/"
    cp -r "$BASE_DIR/resources" "$WORK_DIR/"
    cp -r "$BASE_DIR/gkmsvm_models" "$WORK_DIR/"
    cp "$BASE_DIR/run.sh" "$WORK_DIR/"
    cp "$BASE_DIR/snp_batches/$CHROM.tsv" "$WORK_DIR/input_snp.tsv"

    echo "Starting $CHROM..."
    cd "$WORK_DIR"
    bash run.sh > "$WORK_DIR/run.log" 2>&1 &
    cd "$BASE_DIR"
done

echo "All chromosomes started in parallel!" 
wait
echo "All done!"