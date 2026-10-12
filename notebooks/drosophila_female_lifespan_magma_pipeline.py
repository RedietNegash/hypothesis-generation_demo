#!/usr/bin/env python3

import argparse
import gzip
import os
from pathlib import Path
import re


BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_SUMSTATS = BASE_DIR / "data" / "gwas" / "tmp" / "lifespan_female_gwas_sumstats.tsv"
DEFAULT_GTF = BASE_DIR / "data" / "genes" / "Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
DEFAULT_BFILE = BASE_DIR / "data" / "gwas" / "tmp" / "merged_qc"
DEFAULT_MAGMA_BIN = BASE_DIR / "tools" / "magma" / "magma"
DEFAULT_OUT_DIR = BASE_DIR / "results" / "magma_female"
CHROM_MAP = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "23"}
GENE_ID_PATTERN = re.compile(r'gene_id "([^"]+)"')


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


def write_gene_locations(gtf_path: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f"{output_path.name}.tmp")
    count = 0
    try:
        opener = gzip.open if gtf_path.suffix == ".gz" else open
        with opener(gtf_path, "rt", encoding="utf-8") as source, temporary_path.open(
            "w", encoding="utf-8"
        ) as destination:
            for line in source:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue
                chrom = fields[0].removeprefix("chr")
                if chrom not in CHROM_MAP:
                    continue
                match = GENE_ID_PATTERN.search(fields[8])
                if match is None:
                    raise ValueError(f"Gene record has no gene_id: {line.rstrip()}")
                start, end = int(fields[3]), int(fields[4])
                if start < 1 or end < start or fields[6] not in {"+", "-"}:
                    raise ValueError(f"Invalid gene coordinates or strand: {line.rstrip()}")
                destination.write(
                    f"{match.group(1)}\t{CHROM_MAP[chrom]}\t{start}\t{end}\t{fields[6]}\n"
                )
                count += 1
        if count == 0:
            raise ValueError(f"No Drosophila gene records found in {gtf_path}")
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)
    print(f"Wrote {count} gene locations: {output_path}")
    return output_path


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    validate_inputs(args)
    write_gene_locations(args.gtf, args.out_dir / "gene_loc.txt")


if __name__ == "__main__":
    main()
