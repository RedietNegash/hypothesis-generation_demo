import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PHENOTYPE = REPO_DIR / "data" / "Highfilletal(2016)805pBDSPRRILs.txt"
DEFAULT_GWAS = REPO_DIR / "data" / "Highfill_803RILs_GWAS_PCA.txt"
DEFAULT_GENOTYPE = REPO_DIR / "data" / "Highfill_803RILs_genotype.txt"
DEFAULT_LIFTOVER_CHAIN = REPO_DIR / "data" / "dm3ToDm6.over.chain.gz"
DEFAULT_OUTPUT_DIR = REPO_DIR / "data" / "highfill_finemap"
CHROMOSOMES = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "23"}


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


def lift_dm3_positions(chromosomes, positions, chain_file: Path) -> list:
    from pyliftover import LiftOver

    liftover = LiftOver(str(chain_file))
    mapped_positions = []
    for chrom, pos in zip(chromosomes, positions):
        mapped = liftover.convert_coordinate(f"chr{chrom}", int(pos) - 1) or []
        same_chromosome = [item for item in mapped if item[0] == f"chr{chrom}"]
        mapped_positions.append(
            int(same_chromosome[0][1]) + 1 if len(same_chromosome) == 1 else pd.NA
        )
    return mapped_positions


def prepare_gwas(source_file: Path, chain_file: Path, output_dir: Path) -> Path:
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
    if gwas["CHR"].isna().any() or not gwas["CHR"].isin(CHROMOSOMES).all():
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

    gwas["POS_dm6"] = lift_dm3_positions(gwas["CHR"], gwas["POS"], chain_file)
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


def prepare_genotypes(
    source_file: Path, chain_file: Path, phenotype_file: Path,
    output_dir: Path, plink_bin: str
) -> Path:
    for label, path in (
        ("Highfill genotype", source_file),
        ("dm3-to-dm6 chain", chain_file),
        ("Prepared phenotype", phenotype_file),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} file not found: {path}")
    plink = shutil.which(plink_bin)
    if plink is None:
        raise FileNotFoundError(f"PLINK executable not found: {plink_bin}")

    phenotype = pd.read_csv(phenotype_file, sep="\t", dtype={"patRIL": "string"})
    if "patRIL" not in phenotype or phenotype["patRIL"].duplicated().any():
        raise ValueError("Prepared phenotype requires unique patRIL identifiers")
    ril_ids = phenotype["patRIL"].dropna()

    raw = pd.read_csv(
        source_file, sep=r"\s+", header=None,
        names=["CHR", "POS", "RIL", "A1", "A2", "count1", "count2"],
        dtype={"CHR": "string", "RIL": "string", "A1": "string", "A2": "string"},
    )
    if raw.empty:
        raise ValueError("Highfill genotype file is empty")
    raw["CHR"] = raw["CHR"].str.removeprefix("chr")
    if raw["CHR"].isna().any() or not raw["CHR"].isin(CHROMOSOMES).all():
        raise ValueError("Highfill genotypes contain an unsupported chromosome")
    if raw[["RIL", "A1", "A2"]].isna().any().any() or (raw["A1"] == raw["A2"]).any():
        raise ValueError("Highfill genotypes contain missing or identical alleles")
    for column in ("POS", "count1", "count2"):
        raw[column] = pd.to_numeric(raw[column], errors="raise")
    if not np.isfinite(raw[["POS", "count1", "count2"]].to_numpy()).all():
        raise ValueError("Highfill genotypes contain nonfinite positions or counts")
    if (raw["POS"] <= 0).any() or (raw["POS"] != np.floor(raw["POS"])).any():
        raise ValueError("Genotype positions must be positive integers")
    if (raw[["count1", "count2"]] < 0).any().any():
        raise ValueError("Genotype allele counts must be nonnegative")

    raw["SNP"] = raw["CHR"] + ":" + raw["POS"].astype(int).astype(str)
    if raw.duplicated(["RIL", "SNP"]).any():
        raise ValueError("Highfill genotypes contain duplicate RIL-SNP observations")
    variants = raw[["SNP", "CHR", "POS", "A1", "A2"]].drop_duplicates()
    if variants["SNP"].duplicated().any():
        raise ValueError("Highfill genotypes disagree on alleles for a SNP")
    variants = variants.copy()
    variants["POS_dm6"] = lift_dm3_positions(
        variants["CHR"], variants["POS"], chain_file
    )
    variants = variants.dropna(subset=["POS_dm6"])
    if variants.empty:
        raise ValueError("No genotype variants could be mapped from dm3 to dm6")
    variants["POS_dm6"] = variants["POS_dm6"].astype(int)
    variants["CHR_NUM"] = variants["CHR"].map(CHROMOSOMES)
    variants = variants.sort_values(["CHR_NUM", "POS_dm6", "SNP"])

    missing_rils = ril_ids[~ril_ids.isin(raw["RIL"])]
    if not missing_rils.empty:
        raise ValueError(f"Genotypes are missing for {len(missing_rils)} phenotype RILs")
    ril_ids = ril_ids.reset_index(drop=True)
    if len(ril_ids) < 2:
        raise ValueError("Fewer than two phenotype RILs have genotype data")
    raw = raw.loc[raw["RIL"].isin(ril_ids) & raw["SNP"].isin(variants["SNP"])].copy()
    total = raw["count1"] + raw["count2"]
    fraction_a1 = raw["count1"] / total.replace(0, np.nan)
    raw["call"] = np.where(
        fraction_a1 >= 0.8, 2, np.where(fraction_a1 <= 0.2, 0, np.nan)
    )
    calls = raw.pivot(index="RIL", columns="SNP", values="call")
    calls = calls.reindex(index=ril_ids, columns=variants["SNP"])

    output_dir.mkdir(parents=True, exist_ok=True)
    final_prefix = output_dir / "highfill_genotypes"
    with tempfile.TemporaryDirectory(prefix="highfill_plink_", dir=output_dir) as temp:
        work_dir = Path(temp)
        ped_prefix = work_dir / "highfill_genotypes"
        raw_prefix = work_dir / "highfill_genotypes_raw"
        aligned_prefix = work_dir / "highfill_genotypes_aligned"
        variants[["CHR_NUM", "SNP", "POS_dm6"]].assign(CM=0)[
            ["CHR_NUM", "SNP", "CM", "POS_dm6"]
        ].to_csv(ped_prefix.with_suffix(".map"), sep="\t", header=False, index=False)
        homozygous_a1 = (variants["A1"] + " " + variants["A1"]).to_numpy()
        homozygous_a2 = (variants["A2"] + " " + variants["A2"]).to_numpy()
        with ped_prefix.with_suffix(".ped").open("w") as destination:
            for ril, values in zip(ril_ids, calls.to_numpy()):
                alleles = np.where(
                    values == 2, homozygous_a1,
                    np.where(values == 0, homozygous_a2, "0 0")
                )
                destination.write(f"0\t{ril}\t0\t0\t0\t-9\t" + "\t".join(alleles) + "\n")
        a1_file = work_dir / "effect_alleles.tsv"
        variants[["SNP", "A1"]].to_csv(a1_file, sep="\t", header=False, index=False)

        for command in (
            [plink, "--file", str(ped_prefix), "--make-bed", "--out", str(raw_prefix)],
            [plink, "--bfile", str(raw_prefix), "--a1-allele", str(a1_file),
             "--make-bed", "--out", str(aligned_prefix)],
        ):
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(
                    f"PLINK failed: {' '.join(command)}\n{result.stdout[-2000:]}\n"
                    f"{result.stderr[-2000:]}"
                )

        bim = pd.read_csv(aligned_prefix.with_suffix(".bim"), sep=r"\s+", header=None)
        observed_a1 = bim.set_index(1)[4].astype(str)
        expected_a1 = variants.set_index("SNP")["A1"].astype(str)
        if set(observed_a1.index) != set(expected_a1.index):
            raise ValueError("PLINK reference does not contain all expected variants")
        if not observed_a1.reindex(expected_a1.index).eq(expected_a1).all():
            raise ValueError("PLINK A1 alleles do not match the Highfill genotype source")

        for extension in (".bed", ".bim", ".fam"):
            shutil.copy2(
                aligned_prefix.with_suffix(extension), final_prefix.with_suffix(extension)
            )

    print(f"Prepared {len(ril_ids)} RILs and {len(variants)} variants: {final_prefix}")
    return final_prefix


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Highfill DSPR Drosophila lifespan fine-mapping pipeline"
    )
    parser.add_argument(
        "--stage", choices=("phenotype", "gwas", "genotype"), default="phenotype"
    )
    parser.add_argument("--phenotype", type=Path, default=DEFAULT_PHENOTYPE)
    parser.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE)
    parser.add_argument("--liftover-chain", type=Path, default=DEFAULT_LIFTOVER_CHAIN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--plink-bin", default="plink")
    args = parser.parse_args()

    if args.stage == "phenotype":
        prepare_phenotype(args.phenotype.expanduser(), args.output_dir.expanduser())
    elif args.stage == "gwas":
        prepare_gwas(
            args.gwas.expanduser(), args.liftover_chain.expanduser(),
            args.output_dir.expanduser()
        )
    elif args.stage == "genotype":
        prepare_genotypes(
            args.genotype.expanduser(), args.liftover_chain.expanduser(),
            args.output_dir.expanduser() / "highfill_phenotype_prepared.tsv",
            args.output_dir.expanduser(), args.plink_bin
        )


if __name__ == "__main__":
    main()
