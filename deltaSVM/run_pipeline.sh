#!/bin/bash
cd /mnt/hdd_1/rediet/deltaSVM

for CHR in "$@"; do
    echo "=== Starting $CHR ==="
    bash /mnt/hdd_1/rediet/deltaSVM/restart_chr.sh $CHR
    echo "=== $CHR FULLY DONE ==="
done

echo "ALL DONE!"