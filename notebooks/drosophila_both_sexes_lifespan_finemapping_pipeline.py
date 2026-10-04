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
DEFAULT_GCTA = REPO_DIR / "tools" / "gcta" / "1.94.1" / "gcta64"
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
    raw["A1"] = raw["A1"].str.upper()
    raw["A2"] = raw["A2"].str.upper()
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


def prepare_cojo_input(
    gwas_file: Path, bfile_prefix: Path, output_dir: Path,
    plink_bin: str, p_threshold: float = 1e-5, min_n: int = 100
) -> Path:
    if not 0 < p_threshold < 1 or min_n < 3:
        raise ValueError("COJO requires 0 < p threshold < 1 and minimum N >= 3")
    if not gwas_file.is_file():
        raise FileNotFoundError(f"Prepared Highfill GWAS not found: {gwas_file}")
    for extension in (".bed", ".bim", ".fam"):
        if not bfile_prefix.with_suffix(extension).is_file():
            raise FileNotFoundError(f"PLINK reference is incomplete: {bfile_prefix}")
    plink = shutil.which(plink_bin)
    if plink is None:
        raise FileNotFoundError(f"PLINK executable not found: {plink_bin}")

    gwas = pd.read_csv(gwas_file, sep="\t", dtype={"SNP": "string", "CHR": "string"})
    required = {"SNP", "CHR", "POS_dm6", "BETA", "se", "P", "N"}
    if not required.issubset(gwas.columns):
        raise ValueError(f"Prepared GWAS is missing columns: {sorted(required - set(gwas.columns))}")
    if gwas["SNP"].isna().any() or gwas["SNP"].duplicated().any():
        raise ValueError("Prepared GWAS requires unique SNP identifiers")
    for column in ("POS_dm6", "BETA", "se", "P", "N"):
        gwas[column] = pd.to_numeric(gwas[column], errors="coerce")
    eligible = gwas.loc[(gwas["P"] <= p_threshold) & (gwas["N"] >= min_n)].copy()
    if eligible.empty:
        raise ValueError("No GWAS SNPs pass the COJO p-value and sample-size thresholds")
    if not np.isfinite(eligible[["POS_dm6", "BETA", "se", "P", "N"]].to_numpy()).all():
        raise ValueError("COJO-eligible SNPs contain nonfinite statistics")
    if (eligible["se"] <= 0).any() or (eligible["N"] != np.floor(eligible["N"])).any():
        raise ValueError("COJO-eligible SNPs require positive SE and integer N")

    bim = pd.read_csv(
        bfile_prefix.with_suffix(".bim"), sep=r"\s+", header=None,
        names=["CHR_REF", "SNP", "CM", "POS_REF", "A1_REF", "A2_REF"],
        dtype={"CHR_REF": "string", "SNP": "string", "A1_REF": "string", "A2_REF": "string"},
    )
    if bim["SNP"].duplicated().any():
        raise ValueError("PLINK reference contains duplicate SNP identifiers")

    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="highfill_freq_", dir=output_dir) as temp:
        freq_prefix = Path(temp) / "highfill_genotypes_freq"
        command = [
            plink, "--bfile", str(bfile_prefix), "--keep-allele-order",
            "--freq", "--out", str(freq_prefix),
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(
                f"PLINK frequency calculation failed\n{result.stdout[-2000:]}\n"
                f"{result.stderr[-2000:]}"
            )
        freq = pd.read_csv(freq_prefix.with_suffix(".frq"), sep=r"\s+")

    if freq["SNP"].duplicated().any():
        raise ValueError("PLINK frequency output contains duplicate SNP identifiers")
    reference = bim.merge(
        freq[["SNP", "A1", "A2", "MAF"]], on="SNP", how="left", validate="one_to_one"
    )
    if not reference["A1_REF"].eq(reference["A1"]).all() or not reference["A2_REF"].eq(
        reference["A2"]
    ).all():
        raise ValueError("PLINK frequencies disagree with the reference A1/A2 alleles")

    selected = eligible.merge(reference, on="SNP", how="inner", validate="one_to_one")
    if selected.empty:
        raise ValueError("No COJO-eligible GWAS SNPs are present in the PLINK reference")
    if not selected["CHR"].map(CHROMOSOMES).eq(selected["CHR_REF"]).all():
        raise ValueError("GWAS and reference chromosomes disagree")
    if not selected["POS_dm6"].eq(selected["POS_REF"]).all():
        raise ValueError("GWAS and reference dm6 positions disagree")
    selected["MAF"] = pd.to_numeric(selected["MAF"], errors="coerce")
    selected = selected.loc[selected["MAF"].between(0, 1, inclusive="neither")].copy()
    if selected.empty:
        raise ValueError("No COJO-eligible SNPs have polymorphic reference genotypes")

    cojo = selected[["SNP", "A1_REF", "A2_REF", "MAF", "BETA", "se", "P", "N"]].copy()
    cojo.columns = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
    cojo["N"] = cojo["N"].astype(int)
    cojo = cojo.sort_values(["p", "SNP"])
    output_file = output_dir / "highfill_cojo_input.txt"
    cojo.to_csv(output_file, sep=" ", index=False)
    print(
        f"Prepared {len(cojo)} COJO SNPs from {len(eligible)} eligible GWAS SNPs: "
        f"{output_file}"
    )
    return output_file


def run_cojo(
    gwas_file: Path, bfile_prefix: Path, output_dir: Path,
    plink_bin: str, gcta_bin: str, p_threshold: float = 1e-5,
    min_n: int = 100, collinear: float = 0.5
) -> Path:
    if not 0 < collinear < 1:
        raise ValueError("COJO collinearity cutoff must be between zero and one")
    gcta = shutil.which(gcta_bin)
    if gcta is None:
        raise FileNotFoundError(f"GCTA executable not found: {gcta_bin}")
    cojo_input = prepare_cojo_input(
        gwas_file, bfile_prefix, output_dir, plink_bin, p_threshold, min_n
    )
    result_dir = output_dir / "cojo"
    result_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="highfill_cojo_", dir=result_dir) as temp:
        prefix = Path(temp) / "highfill_lifespan_cojo"
        command = [
            gcta, "--bfile", str(bfile_prefix), "--cojo-file", str(cojo_input),
            "--cojo-slct", "--cojo-p", str(p_threshold),
            "--cojo-collinear", str(collinear), "--out", str(prefix),
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(
                f"GCTA-COJO failed\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
            )
        jma_file = prefix.with_suffix(".jma.cojo")
        if not jma_file.is_file():
            raise RuntimeError("GCTA-COJO did not write a joint-signal result")
        signals = pd.read_csv(jma_file, sep=r"\s+")
        if signals.empty or not {"SNP", "Chr", "bp", "pJ"}.issubset(signals.columns):
            raise ValueError("GCTA-COJO did not return valid independent signals")
        input_snps = set(pd.read_csv(cojo_input, sep=r"\s+")["SNP"])
        if not set(signals["SNP"]).issubset(input_snps):
            raise ValueError("GCTA-COJO selected SNPs absent from its input")
        for artifact in prefix.parent.glob(f"{prefix.name}.*"):
            shutil.copy2(artifact, result_dir / artifact.name)

    output_file = result_dir / "highfill_lifespan_cojo.jma.cojo"
    print(f"COJO selected {len(signals)} independent signals: {output_file}")
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Highfill DSPR Drosophila lifespan fine-mapping pipeline"
    )
    parser.add_argument(
        "--stage", choices=("phenotype", "gwas", "genotype", "cojo"),
        default="phenotype"
    )
    parser.add_argument("--phenotype", type=Path, default=DEFAULT_PHENOTYPE)
    parser.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE)
    parser.add_argument("--liftover-chain", type=Path, default=DEFAULT_LIFTOVER_CHAIN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--plink-bin", default="plink")
    parser.add_argument("--gcta-bin", default=str(DEFAULT_GCTA))
    parser.add_argument("--cojo-p", type=float, default=1e-5)
    parser.add_argument("--min-n", type=int, default=100)
    parser.add_argument("--cojo-collinear", type=float, default=0.5)
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
    elif args.stage == "cojo":
        output_dir = args.output_dir.expanduser()
        run_cojo(
            output_dir / "highfill_gwas_prepared.tsv",
            output_dir / "highfill_genotypes", output_dir,
            args.plink_bin, args.gcta_bin,
            args.cojo_p, args.min_n, args.cojo_collinear,
        )


if __name__ == "__main__":
    main()
