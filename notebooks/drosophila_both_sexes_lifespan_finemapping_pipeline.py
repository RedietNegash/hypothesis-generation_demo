import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
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
DEFAULT_GTF = REPO_DIR / "data" / "genes" / "Drosophila_melanogaster.BDGP6.54.62.chr.gtf.gz"
DEFAULT_OUTPUT_DIR = REPO_DIR / "data" / "highfill_finemap"
DEFAULT_GCTA = REPO_DIR / "tools" / "gcta" / "1.94.1" / "gcta64"
CHROMOSOMES = {"2L": "1", "2R": "2", "3L": "3", "3R": "4", "4": "5", "X": "23"}

QTL_R_CODE = r"""
args <- commandArgs(trailingOnly = TRUE)
phenotype_file <- args[1]
scan_file <- args[2]
peaks_file <- args[3]
raw_file <- args[4]
threshold <- as.numeric(args[5])
reuse_scan <- args[6] == "TRUE"
suppressPackageStartupMessages(library(DSPRqtl))
suppressPackageStartupMessages(library(DSPRqtlDataB))
phenotype <- read.table(phenotype_file, header = TRUE, sep = "\t")
required <- c("patRIL", "Block", "MedLifespanHrs")
if (!all(required %in% names(phenotype))) stop("Prepared phenotype columns are missing")
phenotype$patRIL <- suppressWarnings(as.numeric(phenotype$patRIL))
if (anyNA(phenotype$patRIL) || any(phenotype$patRIL <= 0)) {
  stop("DSPRqtl requires positive numeric patRIL identifiers")
}
if (anyDuplicated(phenotype$patRIL)) stop("Duplicate patRIL identifiers")
if (reuse_scan) {
  scan.results <- readRDS(scan_file)
} else {
  scan.results <- DSPRscan(
    MedLifespanHrs ~ factor(Block), design = "inbredB",
    phenotype.dat = phenotype, id.col = "patRIL"
  )
  temporary_scan <- paste0(scan_file, ".tmp")
  saveRDS(scan.results, temporary_scan)
  if (!file.rename(temporary_scan, scan_file)) stop("Could not save DSPR scan")
}
peaks <- DSPRpeaks(
  scan.results, method = "both", threshold = threshold,
  LODdrop = 2, BCIprob = 0.95
)
saveRDS(peaks, peaks_file)
empty <- data.frame(
  chrom = character(), peak_dm3 = integer(), lod = double(),
  bci_start_dm3 = integer(), bci_end_dm3 = integer(),
  percent_variance = double(), entropy = double(),
  all_founders_observed = logical()
)
if (length(peaks) == 0) {
  raw_peaks <- empty
} else {
  rows <- lapply(peaks, function(item) {
    interval <- item$CI$BCI
    founder_counts <- item$founderNs[paste0("B", 1:8)]
    founder_ok <- length(founder_counts) == 8 &&
      all(is.finite(founder_counts)) && all(founder_counts > 0) &&
      all(is.finite(as.matrix(item$geno.means)))
    data.frame(
      chrom = as.character(item$peak[["chr"]]),
      peak_dm3 = as.integer(item$peak[["Ppos"]]),
      lod = as.numeric(item$peak[["LOD"]]),
      bci_start_dm3 = as.integer(interval[1, "Ppos"]),
      bci_end_dm3 = as.integer(interval[2, "Ppos"]),
      percent_variance = as.numeric(item$perct.var),
      entropy = as.numeric(item$entropy),
      all_founders_observed = founder_ok
    )
  })
  raw_peaks <- do.call(rbind, rows)
}
write.table(raw_peaks, raw_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Found %d raw DSPR peaks at LOD threshold %.2f\n", nrow(raw_peaks), threshold))
"""


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


def collapse_qtl_peaks(peaks: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "chrom", "bci_start_dm3", "bci_end_dm3", "lead_peak_dm3",
        "lead_lod", "raw_peak_count", "all_founders_observed",
    ]
    if peaks.empty:
        return pd.DataFrame(columns=columns)
    required = {
        "chrom", "peak_dm3", "lod", "bci_start_dm3", "bci_end_dm3",
        "all_founders_observed",
    }
    if not required.issubset(peaks.columns):
        raise ValueError(f"DSPR peaks are missing columns: {sorted(required - set(peaks.columns))}")
    peaks = peaks.copy()
    for column in ("peak_dm3", "lod", "bci_start_dm3", "bci_end_dm3"):
        peaks[column] = pd.to_numeric(peaks[column], errors="raise")
    if not np.isfinite(
        peaks[["peak_dm3", "lod", "bci_start_dm3", "bci_end_dm3"]].to_numpy()
    ).all():
        raise ValueError("DSPR peak coordinates and LOD scores must be finite")
    if not (
        (peaks["bci_start_dm3"] <= peaks["peak_dm3"])
        & (peaks["peak_dm3"] <= peaks["bci_end_dm3"])
    ).all():
        raise ValueError("A DSPR peak lies outside its Bayesian interval")

    peaks = peaks.sort_values(["chrom", "bci_start_dm3", "bci_end_dm3", "lod"])
    regions = []
    for row in peaks.itertuples(index=False):
        if (
            regions and regions[-1]["chrom"] == row.chrom
            and row.bci_start_dm3 <= regions[-1]["bci_end_dm3"]
        ):
            region = regions[-1]
            region["bci_end_dm3"] = max(region["bci_end_dm3"], row.bci_end_dm3)
            region["raw_peak_count"] += 1
            region["all_founders_observed"] &= bool(row.all_founders_observed)
            if row.lod > region["lead_lod"]:
                region["lead_peak_dm3"] = row.peak_dm3
                region["lead_lod"] = row.lod
        else:
            regions.append({
                "chrom": row.chrom,
                "bci_start_dm3": row.bci_start_dm3,
                "bci_end_dm3": row.bci_end_dm3,
                "lead_peak_dm3": row.peak_dm3,
                "lead_lod": row.lod,
                "raw_peak_count": 1,
                "all_founders_observed": bool(row.all_founders_observed),
            })
    return pd.DataFrame(regions, columns=columns)


def run_qtl_scan(
    phenotype_file: Path, output_dir: Path, rscript_bin: str,
    threshold: float = 6.8, force_scan: bool = False
) -> Path:
    if not phenotype_file.is_file():
        raise FileNotFoundError(f"Prepared phenotype not found: {phenotype_file}")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("DSPR QTL threshold must be positive and finite")
    rscript = shutil.which(rscript_bin)
    if rscript is None:
        raise FileNotFoundError(f"Rscript executable not found: {rscript_bin}")

    qtl_dir = output_dir / "dspr_qtl"
    qtl_dir.mkdir(parents=True, exist_ok=True)
    scan_file = qtl_dir / "dspr_scan.rds"
    peaks_file = qtl_dir / "dspr_peaks.rds"
    raw_file = qtl_dir / "dspr_raw_peaks.tsv"
    regions_file = qtl_dir / "dspr_qtl_regions.tsv"
    manifest_file = qtl_dir / "scan_manifest.json"
    phenotype_hash = hashlib.sha256(phenotype_file.read_bytes()).hexdigest()
    manifest = {
        "phenotype_sha256": phenotype_hash,
        "scan_model": "DSPRscan_inbredB_MedLifespanHrs_factorBlock_v1",
    }
    try:
        cached_manifest = json.loads(manifest_file.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        cached_manifest = None
    reuse_scan = not force_scan and scan_file.is_file() and cached_manifest == manifest

    command = [
        rscript, "--vanilla", "-e", QTL_R_CODE,
        str(phenotype_file), str(scan_file), str(peaks_file), str(raw_file),
        str(threshold), "TRUE" if reuse_scan else "FALSE",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"DSPR QTL scan failed\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
        )
    if not raw_file.is_file() or not peaks_file.is_file() or not scan_file.is_file():
        raise RuntimeError("DSPR QTL scan did not create its expected outputs")
    peaks = pd.read_csv(raw_file, sep="\t")
    regions = collapse_qtl_peaks(peaks)
    regions.to_csv(regions_file, sep="\t", index=False)
    manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
    print(result.stdout.strip())
    print(f"Collapsed {len(peaks)} raw peaks into {len(regions)} QTL regions: {regions_file}")
    return regions_file


def load_gene_annotations(gtf_file: Path) -> pd.DataFrame:
    if not gtf_file.is_file():
        raise FileNotFoundError(f"Drosophila gene annotation not found: {gtf_file}")
    rows = []
    with gzip.open(gtf_file, "rt", encoding="utf-8") as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            gene_id = attributes.get("gene_id", "")
            rows.append({
                "chrom": fields[0].removeprefix("chr"),
                "start": int(fields[3]),
                "end": int(fields[4]),
                "gene_id": gene_id,
                "gene_name": attributes.get("gene_name", gene_id),
                "gene_biotype": attributes.get("gene_biotype", ""),
            })
    genes = pd.DataFrame(rows)
    if genes.empty or genes["gene_id"].eq("").any():
        raise ValueError(f"No valid gene features found in {gtf_file}")
    return genes


def map_qtl_peaks_to_genes(
    regions_file: Path, chain_file: Path, gtf_file: Path, output_dir: Path
) -> Path:
    from pyliftover import LiftOver

    if not regions_file.is_file():
        raise FileNotFoundError(f"DSPR QTL regions not found: {regions_file}")
    if not chain_file.is_file():
        raise FileNotFoundError(f"dm3-to-dm6 chain file not found: {chain_file}")
    regions = pd.read_csv(regions_file, sep="\t", dtype={"chrom": "string"})
    required = {
        "chrom", "bci_start_dm3", "bci_end_dm3", "lead_peak_dm3",
        "lead_lod", "raw_peak_count", "all_founders_observed",
    }
    if not required.issubset(regions.columns):
        raise ValueError(f"DSPR regions are missing columns: {sorted(required - set(regions.columns))}")
    genes = load_gene_annotations(gtf_file)
    liftover = LiftOver(str(chain_file))
    columns = [
        "chrom", "bci_start_dm3", "bci_end_dm3", "lead_peak_dm3",
        "lead_peak_dm6", "lead_lod", "raw_peak_count", "all_founders_observed",
        "gene_id", "gene_name", "gene_biotype", "relation", "distance_bp",
    ]
    rows = []
    for region in regions.itertuples(index=False):
        chrom = str(region.chrom).removeprefix("chr")
        peak_dm3 = int(region.lead_peak_dm3)
        if chrom not in CHROMOSOMES or peak_dm3 < 1:
            raise ValueError(f"Invalid QTL peak coordinate: {chrom}:{peak_dm3}")
        mapped = liftover.convert_coordinate(f"chr{chrom}", peak_dm3 - 1) or []
        mapped = [hit for hit in mapped if hit[0] == f"chr{chrom}"]
        if len(mapped) != 1:
            raise ValueError(f"QTL peak has no unique same-chromosome dm6 mapping: {chrom}:{peak_dm3}")
        peak_dm6 = int(mapped[0][1]) + 1
        on_chrom = genes.loc[genes["chrom"] == chrom]
        if on_chrom.empty:
            raise ValueError(f"No dm6 genes found on chromosome {chrom}")
        overlapping = on_chrom.loc[
            (on_chrom["start"] <= peak_dm6) & (peak_dm6 <= on_chrom["end"])
        ]
        if overlapping.empty:
            distance = np.maximum(
                peak_dm6 - on_chrom["end"], on_chrom["start"] - peak_dm6
            )
            nearest_distance = int(distance.min())
            selected = on_chrom.loc[distance == nearest_distance]
            relation = "nearest"
        else:
            selected = overlapping
            nearest_distance = 0
            relation = "overlapping"
        for gene in selected.itertuples(index=False):
            rows.append({
                "chrom": chrom,
                "bci_start_dm3": int(region.bci_start_dm3),
                "bci_end_dm3": int(region.bci_end_dm3),
                "lead_peak_dm3": peak_dm3,
                "lead_peak_dm6": peak_dm6,
                "lead_lod": float(region.lead_lod),
                "raw_peak_count": int(region.raw_peak_count),
                "all_founders_observed": bool(region.all_founders_observed),
                "gene_id": gene.gene_id,
                "gene_name": gene.gene_name,
                "gene_biotype": gene.gene_biotype,
                "relation": relation,
                "distance_bp": nearest_distance,
            })
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "highfill_qtl_peak_genes.tsv"
    pd.DataFrame(rows, columns=columns).to_csv(output_file, sep="\t", index=False)
    print(f"Mapped {len(regions)} DSPR QTL peaks to {len(rows)} gene records: {output_file}")
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Highfill DSPR Drosophila lifespan fine-mapping pipeline"
    )
    parser.add_argument(
        "--stage", choices=("phenotype", "gwas", "genotype", "cojo", "qtl", "genes"),
        default="phenotype"
    )
    parser.add_argument("--phenotype", type=Path, default=DEFAULT_PHENOTYPE)
    parser.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE)
    parser.add_argument("--liftover-chain", type=Path, default=DEFAULT_LIFTOVER_CHAIN)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--plink-bin", default="plink")
    parser.add_argument("--gcta-bin", default=str(DEFAULT_GCTA))
    parser.add_argument("--cojo-p", type=float, default=1e-5)
    parser.add_argument("--min-n", type=int, default=100)
    parser.add_argument("--cojo-collinear", type=float, default=0.5)
    parser.add_argument("--rscript-bin", default="Rscript")
    parser.add_argument("--qtl-threshold", type=float, default=6.8)
    parser.add_argument("--force-qtl-scan", action="store_true")
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
    elif args.stage == "qtl":
        output_dir = args.output_dir.expanduser()
        run_qtl_scan(
            output_dir / "highfill_phenotype_prepared.tsv", output_dir,
            args.rscript_bin, args.qtl_threshold, args.force_qtl_scan,
        )
    elif args.stage == "genes":
        output_dir = args.output_dir.expanduser()
        map_qtl_peaks_to_genes(
            output_dir / "dspr_qtl" / "dspr_qtl_regions.tsv",
            args.liftover_chain.expanduser(), args.gtf.expanduser(), output_dir,
        )


if __name__ == "__main__":
    main()
