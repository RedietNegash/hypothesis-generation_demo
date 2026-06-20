import pandas as pd
import glob
import os

print("Loading OpenTargets variants...")
files = glob.glob('/mnt/hdd_1/rediet/opentargets_variants/variant/part-*.parquet')
ot_variants = set()
for f in files:
    df = pd.read_parquet(f, columns=['variantId'])
    df['variantId'] = 'chr' + df['variantId']
    ot_variants.update(df['variantId'].tolist())
    print(f"Loaded {f.split('/')[-1]}: {len(ot_variants):,} variants so far")

print(f"\nTotal OT variants: {len(ot_variants):,}")

for CHR in ['chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']:
    snp_file = f'/mnt/hdd_1/rediet/deltaSVM/snp_batches/{CHR}.tsv'
    out_file = f'/mnt/hdd_1/rediet/deltaSVM/snp_batches/{CHR}_ot.tsv'

    with open(snp_file) as f:
        snps = [line.strip() for line in f]

    filtered = [s for s in snps if s in ot_variants]

    with open(out_file, 'w') as f:
        f.write('\n'.join(filtered))

    print(f"{CHR}: {len(snps):,} → {len(filtered):,} SNPs after filtering")
