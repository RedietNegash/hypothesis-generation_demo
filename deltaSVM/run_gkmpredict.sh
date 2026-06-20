#!/bin/bash

BASE_DIR="/mnt/hdd_1/rediet/deltaSVM"

run_gkm_for_chr() {
    CHR=$1
    WORK_DIR="$BASE_DIR/runs/$CHR"
    cd "$WORK_DIR"

    echo "[$CHR] Starting gkmpredict..."
    while [ ! -s out/pbs.pred.tsv ]; do
        sleep 60
    done

    sort -k1,1 -k3,3 -o out/pbs.pred.tsv out/pbs.pred.tsv

    cat resources/thresholds.obs.tsv | awk '{print $1" "$3}' | while read tf model; do
        if [ ! -f "tmp/${tf}.merge.gkm.tsv" ]; then echo "$tf $model"; fi
    done > missing_tfs.txt

    echo "[$CHR] $(wc -l < missing_tfs.txt) TFs to process"

    cat missing_tfs.txt | parallel -j4 --colsep ' ' '
        tf={1}; model={2}
        scripts/gkmpredict -T 16 data/selex_allelic_oligos.ref.fa gkmsvm_models/${model}.model.txt tmp/${tf}.ref.gkm.tsv
        sort -k1,1 -o tmp/${tf}.ref.gkm.tsv tmp/${tf}.ref.gkm.tsv
        scripts/gkmpredict -T 16 data/selex_allelic_oligos.alt.fa gkmsvm_models/${model}.model.txt tmp/${tf}.alt.gkm.tsv
        sort -k1,1 -o tmp/${tf}.alt.gkm.tsv tmp/${tf}.alt.gkm.tsv
        paste tmp/${tf}.ref.gkm.tsv tmp/${tf}.alt.gkm.tsv | cut -f1,2,4 > tmp/${tf}.merge.gkm.tsv
        echo "Done: $tf"
    '

    echo "[$CHR] gkmpredict done!"
}

export -f run_gkm_for_chr
export BASE_DIR

for CHR in chr15 chr16 chr17 chr18 chr19 chr20; do
    run_gkm_for_chr $CHR &
done

wait
echo "All chromosomes done!"