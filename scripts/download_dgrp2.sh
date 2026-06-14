#!/usr/bin/env bash
set -euo pipefail

OUTDIR="data/reference"
TMP="data/reference/tmp"
CHROMS=(2L 2R 3L 3R 4 X)

BASE_URL="https://zenodo.org/records/837947/files"
STEM="dgrp2_dm6_dbSNP.vcf"

mkdir -p "$OUTDIR" "$TMP"

for EXT in bed bim fam; do
    DEST="$TMP/${STEM}.${EXT}"
    if [[ -f "$DEST" ]]; then
        echo "${EXT}: already downloaded"
    else
        echo "Downloading ${STEM}.${EXT} ..."
        wget -c -O "$DEST" "${BASE_URL}/${STEM}.${EXT}?download=1"
    fi
done

echo ""
echo "Chromosome arms found in bim:"
awk '{print $1}' "$TMP/${STEM}.bim" | sort -u

echo ""
for CHROM in "${CHROMS[@]}"; do
    OUT="$OUTDIR/DGRP.$CHROM"
    if [[ -f "${OUT}.bed" && -f "${OUT}.bim" && -f "${OUT}.fam" ]]; then
        echo "chr$CHROM: plink files already exist, skipping"
        continue
    fi
    echo "chr$CHROM: extracting ..."
    plink --bfile "$TMP/$STEM" \
          --chr "$CHROM" \
          --allow-extra-chr \
          --keep-allele-order \
          --make-bed \
          --out "$OUT" \
          --silent
    N=$(wc -l < "${OUT}.bim")
    echo "  -> ${N} SNPs written"
done

echo ""
echo "Done. Plink files in $OUTDIR/"
ls -lh "$OUTDIR"/*.bim 2>/dev/null || true

read -rp "Remove temp merged plink files to save space? [y/N] " RESP
if [[ "${RESP,,}" == "y" ]]; then
    rm -rf "$TMP"
    echo "Temp files removed."
fi
