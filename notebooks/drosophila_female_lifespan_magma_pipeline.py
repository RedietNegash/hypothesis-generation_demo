#!/usr/bin/env python3

import argparse
import gzip
import os
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_SUMSTATS = BASE_DIR / "data" / "gwas" / "tmp" / "lifespan_female_gwas_sumstats.tsv"
DEFAULT_GTF = BASE_DIR / "data" / "genes" / "Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
DEFAULT_BFILE = BASE_DIR / "data" / "gwas" / "tmp" / "merged_qc"
DEFAULT_MAGMA_BIN = BASE_DIR / "tools" / "magma" / "magma"
DEFAULT_OUT_DIR = BASE_DIR / "results" / "magma_female"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MAGMA gene and gene-set analysis for DGRP female lifespan."
    )
    parser.add_argument("--sumstats", type=Path, default=DEFAULT_SUMSTATS)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--bfile", type=Path, default=DEFAULT_BFILE)
    parser.add_argument("--magma-bin", type=Path, default=DEFAULT_MAGMA_BIN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args(argv)


def validate_inputs(args: argparse.Namespace) -> None:
    required_files = [
        args.sumstats,
        args.gtf,
        args.magma_bin,
        *(Path(f"{args.bfile}{extension}") for extension in (".bed", ".bim", ".fam")),
    ]
    missing = [path for path in required_files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing MAGMA inputs:\n" + "\n".join(str(path) for path in missing)
        )
    if not os.access(args.magma_bin, os.X_OK):
        raise PermissionError(f"MAGMA binary is not executable: {args.magma_bin}")

    opener = gzip.open if args.sumstats.suffix == ".gz" else open
    with opener(args.sumstats, "rt", encoding="utf-8") as source:
        columns = source.readline().split()
        required_columns = {"CHR", "POS", "SNP", "p", "N"}
        missing_columns = required_columns - set(columns)
        if missing_columns:
            raise ValueError(
                f"GWAS summary statistics are missing columns: {sorted(missing_columns)}"
            )
        if not any(line.strip() for line in source):
            raise ValueError(f"GWAS summary statistics contain no variants: {args.sumstats}")

    print(f"Validated MAGMA inputs: {args.sumstats}")


def main(argv: list[str] | None = None) -> None:
    validate_inputs(parse_args(argv))


if __name__ == "__main__":
    main()
