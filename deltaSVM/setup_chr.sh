#!/bin/bash
CHR=$1
BASE_DIR="/mnt/hdd_1/rediet/deltaSVM"

echo "Setting up $CHR..."
mkdir -p $BASE_DIR/runs/$CHR/{data,tmp,out,log}
cp -r $BASE_DIR/scripts $BASE_DIR/runs/$CHR/
cp -r $BASE_DIR/resources $BASE_DIR/runs/$CHR/
cp -r $BASE_DIR/gkmsvm_models $BASE_DIR/runs/$CHR/
cp $BASE_DIR/snp_batches/$CHR.tsv $BASE_DIR/runs/$CHR/input_snp.tsv

cd $BASE_DIR/runs/$CHR
python scripts/generate_allelic_seqs.py -f resources/hs38/hs38.fa -s input_snp.tsv -o data/selex_allelic_oligos
scripts/deltasvm_subset_multi data/selex_allelic_oligos.ref.fa data/selex_allelic_oligos.alt.fa resources/models.weights.txt out/pbs.pred.tsv resources/thresholds.pbs.tsv
echo "$CHR setup done!"