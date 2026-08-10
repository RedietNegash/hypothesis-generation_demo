#!/usr/bin/env python3
"""
verify_harmonisation.py

Verify that Open Targets / EBI harmonised GWAS files can be used directly
in the LDSC h2-cts pipeline without an additional harmonisation step.

Steps:
  1. Sample studies from Open Targets parquet (FinnGen + GCST with QC values)
  2. Download harmonised files from EBI FTP
  3. Run LDSC with harmonised columns (hm_*) vs raw columns
  4. Compare results
"""
