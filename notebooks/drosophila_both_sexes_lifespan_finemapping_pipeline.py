import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PHENOTYPE = REPO_DIR / "data" / "Highfilletal(2016)805pBDSPRRILs.txt"
DEFAULT_GWAS = REPO_DIR / "data" / "Highfill_803RILs_GWAS_PCA.txt"
DEFAULT_LIFTOVER_CHAIN = REPO_DIR / "data" / "dm3ToDm6.over.chain.gz"
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


def prepare_gwas(source_file: Path, chain_file: Path, output_dir: Path) -> Path:
    from pyliftover import LiftOver
    from scipy.stats import t as t_distribution

    if not source_file.is_file():
        raise FileNotFoundError(f"Highfill GWAS file not found: {source_file}")
    if not chain_file.is_file():
        raise FileNotFoundError(f"dm3-to-dm6 chain file not found: {chain_file}")

    gwas = pd.read_csv(source_file, sep=r"\s+", dtype={"SNP": "string", "CHR": "string"})
    required = {"SNP", "CHR", "POS", "BETA", "P", "N"}
    missing = required - set(gwas.columns)
    if missing:
        raise ValueError(f"Highfill GWAS is missing columns: {sorted(missing)}")
    if gwas.empty or gwas["SNP"].isna().any() or gwas["SNP"].duplicated().any():
        raise ValueError("Highfill GWAS requires unique, nonmissing SNP identifiers")

    gwas["CHR"] = gwas["CHR"].str.removeprefix("chr")
    if gwas["CHR"].isna().any() or not gwas["CHR"].isin(
        ["2L", "2R", "3L", "3R", "4", "X"]
    ).all():
        raise ValueError("Highfill GWAS contains an unsupported chromosome")
    for column in ("POS", "BETA", "P", "N"):
        gwas[column] = pd.to_numeric(gwas[column], errors="raise")
    if not np.isfinite(gwas[["POS", "BETA", "P", "N"]].to_numpy()).all():
        raise ValueError("Highfill GWAS contains nonfinite statistics")
    if (gwas["POS"] <= 0).any() or (gwas["POS"] != np.floor(gwas["POS"])).any():
        raise ValueError("GWAS positions must be positive integers")
    if (gwas["N"] <= 2).any() or (gwas["N"] != np.floor(gwas["N"])).any():
        raise ValueError("GWAS sample sizes must be integers greater than two")
    if ((gwas["P"] < 0) | (gwas["P"] > 1)).any():
        raise ValueError("GWAS p-values must be between zero and one")

    liftover = LiftOver(str(chain_file))
    positions = []
    for chrom, pos in zip(gwas["CHR"], gwas["POS"]):
        mapped = liftover.convert_coordinate(f"chr{chrom}", int(pos) - 1) or []
        same_chromosome = [item for item in mapped if item[0] == f"chr{chrom}"]
        positions.append(int(same_chromosome[0][1]) + 1 if len(same_chromosome) == 1 else pd.NA)
    gwas["POS_dm6"] = positions
    unmapped = int(gwas["POS_dm6"].isna().sum())
    gwas = gwas.dropna(subset=["POS_dm6"]).copy()
    if gwas.empty:
        raise ValueError("No GWAS variants could be mapped from dm3 to dm6")
    gwas["POS_dm6"] = gwas["POS_dm6"].astype(int)

    p_values = gwas["P"].clip(lower=1e-300)
    t_statistic = t_distribution.isf(p_values / 2, df=gwas["N"] - 2)
    valid_se = (gwas["BETA"] != 0) & np.isfinite(t_statistic) & (t_statistic > 0)
    gwas["se"] = np.nan
    gwas.loc[valid_se, "se"] = np.abs(gwas.loc[valid_se, "BETA"]) / t_statistic[valid_se]

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "highfill_gwas_prepared.tsv"
    gwas.to_csv(output_file, sep="\t", index=False)
    print(f"Prepared {len(gwas)} GWAS variants; {unmapped} unmapped: {output_file}")
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Highfill DSPR Drosophila lifespan fine-mapping pipeline"
    )
    parser.add_argument("--stage", choices=("phenotype", "gwas"), default="phenotype")
    parser.add_argument("--phenotype", type=Path, default=DEFAULT_PHENOTYPE)
    parser.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    parser.add_argument("--liftover-chain", type=Path, default=DEFAULT_LIFTOVER_CHAIN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    if args.stage == "phenotype":
        prepare_phenotype(args.phenotype.expanduser(), args.output_dir.expanduser())
    elif args.stage == "gwas":
        prepare_gwas(
            args.gwas.expanduser(), args.liftover_chain.expanduser(),
            args.output_dir.expanduser()
        )


if __name__ == "__main__":
    main()
