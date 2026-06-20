import gzip
import sys

input_vcf = "/mnt/hdd_1/abdu_md/biocypher_data_bizon/dbsnp/00-common_all.vcf.gz"
output_snp = "input_snp_common.tsv"
output_rsid_map = "rsid_map.tsv"

count = 0
written = 0

with gzip.open(input_vcf, "rt") as f, \
     open(output_snp, "w") as out_snp, \
     open(output_rsid_map, "w") as out_map:

    out_map.write("variant\trsId\n")

    for line in f:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")e 
        chrom, pos, rsid, ref, alt, info = cols[0], cols[1], cols[2], cols[3], cols[4], cols[7]

        count += 1
        if "VC=SNV" not in info:
            continue
        if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
            continue

       
        if len(ref) != 1 or len(alt) != 1:
            continue

        variant = f"chr{chrom}_{pos}_{ref}_{alt}"
        out_snp.write(variant + "\n")
        out_map.write(f"{variant}\t{rsid}\n")
        written += 1

        if written % 1000000 == 0:
            print(f"Processed {count:,} lines, written {written:,} SNPs...")

print(f"Done! Total lines: {count:,}, SNPs written: {written:,}")