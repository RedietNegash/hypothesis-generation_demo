#!/bin/bash
CHR=$1
cd /mnt/hdd_1/rediet/deltaSVM/runs/$CHR

cat resources/thresholds.obs.tsv | awk '{print $1" "$3}' | while read tf model; do
    if [ ! -f "tmp/${tf}.merge.gkm.tsv" ]; then echo "$tf $model"; fi
done > missing_tfs.txt
echo "$CHR missing: $(wc -l < missing_tfs.txt) TFs"

cat missing_tfs.txt | parallel -j4 --colsep ' ' '
    tf={1}; model={2}
    scripts/gkmpredict -T 16 data/selex_allelic_oligos.ref.fa gkmsvm_models/${model}.model.txt tmp/${tf}.ref.gkm.tsv
    sort -k1,1 -o tmp/${tf}.ref.gkm.tsv tmp/${tf}.ref.gkm.tsv
    scripts/gkmpredict -T 16 data/selex_allelic_oligos.alt.fa gkmsvm_models/${model}.model.txt tmp/${tf}.alt.gkm.tsv
    sort -k1,1 -o tmp/${tf}.alt.gkm.tsv tmp/${tf}.alt.gkm.tsv
    paste tmp/${tf}.ref.gkm.tsv tmp/${tf}.alt.gkm.tsv | cut -f1,2,4 > tmp/${tf}.merge.gkm.tsv
    echo "Done: $tf"
'

python scripts/obs_pred.py -s resources/thresholds.obs.tsv -o out/obs.pred.tsv  
sort -k1,1 -k2,2 -o out/obs.pred.tsv out/obs.pred.tsv
paste out/obs.pred.tsv out/pbs.pred.tsv | cut -f1-5,7,9 | sort -k7,7 -k5,5r | sed '1i snp\ttf\tallele1_bind\tallele2_bind\tseq_binding\tdeltaSVM_score\tpreferred_allele' > out/summary.pred.tsv

python3 -c "
import pandas as pd
df = pd.read_csv('out/summary.pred.tsv', sep='\t', low_memory=False)
df = df.rename(columns={'snp':'variant','tf':'TfName','preferred_allele':'Effect','deltaSVM_score':'Score'})
df['rsId'] = '.'
df['TfId'] = '.'
df = df[df['Effect'].isin(['Gain','Loss'])]
df = df[['variant','rsId','TfId','TfName','Effect','Score']]
df.to_parquet('out/${CHR}_deltasvm.parquet', index=False)
print('${CHR} parquet done! rows:', len(df))
"
echo "=== $CHR FULLY DONE ==="