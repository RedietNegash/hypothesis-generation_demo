import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PHENOTYPE = REPO_DIR / "data" / "Highfilletal(2016)805pBDSPRRILs.txt"
DEFAULT_OUTPUT_DIR = REPO_DIR / "data" / "highfill_finemap"


def prepare_phenotype(source_file: Path, output_dir: Path) -> Path:
    if not source_file.is_file():
        raise FileNotFoundError(f"Highfill phenotype file not found: {source_file}")

    source = pd.read_csv(source_file, sep=r"\s+", dtype={"RIL": "string", "Block": "string"})
    required = {"RIL", "Block", "MedLifespanHrs"}
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"Highfill phenotype is missing columns: {sorted(missing)}")

    phenotype = source[["RIL", "Block", "MedLifespanHrs"]].copy()
    phenotype["RIL"] = phenotype["RIL"].str.strip()
    phenotype["Block"] = phenotype["Block"].str.strip()
    if phenotype["RIL"].isna().any() or phenotype["RIL"].eq("").any():
        raise ValueError("Highfill phenotype contains missing RIL identifiers")
    if phenotype["RIL"].duplicated().any():
        raise ValueError("Highfill phenotype contains duplicate RIL identifiers")
    if phenotype["Block"].isna().any() or phenotype["Block"].eq("").any():
        raise ValueError("Highfill phenotype contains missing block assignments")

    phenotype["MedLifespanHrs"] = pd.to_numeric(
        phenotype["MedLifespanHrs"], errors="raise"
    )
    phenotype = phenotype.dropna(subset=["MedLifespanHrs"])
    if len(phenotype) < 2 or not np.isfinite(phenotype["MedLifespanHrs"]).all():
        raise ValueError("At least two finite lifespan values are required")
    if (phenotype["MedLifespanHrs"] <= 0).any():
        raise ValueError("Lifespan values must be positive")

    phenotype = phenotype.rename(columns={"RIL": "patRIL"})
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "highfill_phenotype_prepared.tsv"
    phenotype.to_csv(output_file, sep="\t", index=False)
    print(f"Prepared {len(phenotype)} of {len(source)} RIL phenotypes: {output_file}")
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Highfill DSPR Drosophila lifespan fine-mapping pipeline"
    )
    parser.add_argument("--stage", choices=("phenotype",), default="phenotype")
    parser.add_argument("--phenotype", type=Path, default=DEFAULT_PHENOTYPE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    if args.stage == "phenotype":
        prepare_phenotype(args.phenotype.expanduser(), args.output_dir.expanduser())


if __name__ == "__main__":
    main()
