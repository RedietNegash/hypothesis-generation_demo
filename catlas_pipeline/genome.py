import os
import subprocess

import pandas as pd
from pyfaidx import Fasta

UCSC_FASTA_URLS = {
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
}


def load_genome(fasta_path: str) -> Fasta:
    return Fasta(fasta_path)


def ensure_reference_genomes(paths) -> None:
    """Make sure hg38.fa and hg19.fa exist under paths.scratch.

    Both builds are needed: hg19 for the reference LD panel / build detection,
    hg38 as the liftover target for sequence extraction. For each, if it is not
    already in scratch we first symlink a copy sitting in the project base
    (``{base}/hg38.fa``), otherwise download+decompress it from UCSC (~1 GB each).
    """
    os.makedirs(paths.scratch, exist_ok=True)
    for build in ("hg38", "hg19"):
        fa = getattr(paths, f"fasta_{build}")
        if os.path.exists(fa):
            continue
        repo_copy = f"{paths.base}/{build}.fa"
        if os.path.exists(repo_copy):
            print(f"  Linking {build}.fa from {repo_copy}")
            os.symlink(repo_copy, fa)
            continue
        gz = f"{fa}.gz"
        print(f"  Downloading {build}.fa from UCSC (~1 GB, one-time)...")
        subprocess.run(
            [
                "wget",
                "--continue",
                "--tries=5",
                "--timeout=120",
                "-O",
                gz,
                UCSC_FASTA_URLS[build],
            ],
            check=True,
        )
        subprocess.run(["gunzip", "-f", gz], check=True)
        print(f"  Ready: {fa}")


def detect_genome_build(
    df: pd.DataFrame,
    genome_hg19: Fasta,
    genome_hg38: Fasta,
    n_test: int = 50,
    seed: int = 0,
) -> str:

    sample = df.sample(min(n_test, len(df)), random_state=seed)

    def _match_rate(genome: Fasta) -> float:

        matches = 0
        checked = 0
        for _, row in sample.iterrows():
            chrom = (
                f"chr{row['chr']}"
                if not str(row["chr"]).startswith("chr")
                else str(row["chr"])
            )
            pos = int(row["pos"])
            try:
                seq = genome[chrom][pos - 1 : pos].seq.upper()
            except Exception:
                continue
            checked += 1
            if seq == row["ref"] or seq == row["alt"]:
                matches += 1
        return matches / checked if checked else 0.0

    rate_hg19 = _match_rate(genome_hg19)
    rate_hg38 = _match_rate(genome_hg38)
    print(
        f"Build detection: hg19 ref/alt-match={rate_hg19:.1%}, "
        f"hg38 ref/alt-match={rate_hg38:.1%} (n={len(sample)} sampled variants)"
    )

    HIGH, LOW, GAP = 0.85, 0.30, 0.30
    if rate_hg19 >= HIGH and rate_hg19 - rate_hg38 >= GAP:
        return "hg19"
    if rate_hg38 >= HIGH and rate_hg38 - rate_hg19 >= GAP:
        return "hg38"
    raise ValueError(
        f"Ambiguous genome build: hg19={rate_hg19:.1%}, hg38={rate_hg38:.1%} match rate. "
        f"Neither build clearly wins -- check ref/alt column mapping and chr/pos parsing "
        f"in this dataset's prep function before proceeding. Do not guess."
    )


def liftover_hg19_to_hg38(
    df_in: pd.DataFrame,
    liftover_bin: str,
    chain_file: str,
    dataset_dir: str,
    chrom_col="chr",
    pos_col="pos",
    rsid_col="rsid",
    tag="variant",
) -> pd.DataFrame:

    df_in = df_in.drop_duplicates(subset=rsid_col, keep="first").copy()

    bed_in = f"{dataset_dir}/{tag}_hg19.bed"
    with open(bed_in, "w") as f:
        for _, row in df_in.iterrows():
            chrom = (
                f"chr{row[chrom_col]}"
                if not str(row[chrom_col]).startswith("chr")
                else str(row[chrom_col])
            )
            pos = int(row[pos_col])
            f.write(f"{chrom}\t{pos - 1}\t{pos}\t{row[rsid_col]}\n")

    bed_out = f"{dataset_dir}/{tag}_hg38.bed"
    unmapped = f"{dataset_dir}/{tag}_unmapped.bed"
    result = subprocess.run(
        [liftover_bin, bed_in, chain_file, bed_out, unmapped],
        capture_output=True,
        text=True,
    )
    print(f"  liftOver stderr: {result.stderr.strip()}")

    lifted = pd.read_csv(
        bed_out, sep="\t", header=None, names=["chr_lifted", "start", "end", rsid_col]
    )
    lifted["pos_hg38"] = lifted["end"].astype("Int64")
    lifted = lifted.drop_duplicates(subset=rsid_col, keep="first")

    print(f"  Before liftover: {len(df_in):,}  |  Successfully lifted: {len(lifted):,}")

    df_out = df_in.drop(columns=[pos_col]).merge(
        lifted[[rsid_col, "pos_hg38"]].rename(columns={"pos_hg38": pos_col}),
        on=rsid_col,
        how="inner",
    )
    df_out[chrom_col] = df_out[chrom_col].astype(str)
    df_out[pos_col] = df_out[pos_col].astype(int)
    print(f"  After merge back: {len(df_out):,} variants retained")
    return df_out


def verify_against_hg38(
    df_in: pd.DataFrame,
    genome_hg38: Fasta,
    n_test: int = 20,
    chrom_col="chr",
    pos_col="pos",
    ref_col="ref",
    alt_col="alt",
    label="",
    seed: int = 0,
) -> int:
    """Spot-check that ref/alt line up with hg38.fa at the given locus (build/liftover
    correctness). Matches against {ref, alt} since some GWAS files label alleles as
    effect/non-effect rather than genomic ref/alt -- allele orientation itself is
    fixed afterwards by reorient_alleles. Returns mismatch count out of n_test."""
    sample = df_in.sample(min(n_test, len(df_in)), random_state=seed)
    mismatches = 0
    for _, row in sample.iterrows():
        chrom = (
            f"chr{row[chrom_col]}"
            if not str(row[chrom_col]).startswith("chr")
            else str(row[chrom_col])
        )
        pos = int(row[pos_col])
        seq = genome_hg38[chrom][pos - 1 : pos].seq.upper()
        if seq != row[ref_col] and seq != row[alt_col]:
            mismatches += 1
    print(f"{label} mismatches: {mismatches}/{len(sample)}")
    return mismatches


def reorient_alleles(
    df_in: pd.DataFrame,
    genome_hg38: Fasta,
    chrom_col="chr",
    pos_col="pos",
    ref_col="ref",
    alt_col="alt",
    beta_col="beta",
) -> pd.DataFrame:

    def _fix(row):
        chrom = (
            f"chr{row[chrom_col]}"
            if not str(row[chrom_col]).startswith("chr")
            else str(row[chrom_col])
        )
        pos = int(row[pos_col])
        seq = genome_hg38[chrom][pos - 1 : pos].seq.upper()
        if seq == row[ref_col]:
            return row[ref_col], row[alt_col], row[beta_col]
        elif seq == row[alt_col]:
            return row[alt_col], row[ref_col], -row[beta_col]
        else:
            return None, None, None

    results = df_in.apply(_fix, axis=1, result_type="expand")
    results.columns = ["ref_fixed", "alt_fixed", "beta_fixed"]
    df_out = df_in.copy()
    df_out[ref_col] = results["ref_fixed"]
    df_out[alt_col] = results["alt_fixed"]
    df_out[beta_col] = results["beta_fixed"]

    before = len(df_out)
    df_out = df_out.dropna(subset=[ref_col, alt_col])
    print(f"  Dropped {before - len(df_out)} true mismatches during reorientation")
    return df_out


def extract_sequences(
    df_variants: pd.DataFrame, genome_hg38: Fasta, half_window: int, out_path: str
) -> pd.DataFrame:
    rows, skipped = [], 0
    for _, row in df_variants.iterrows():
        chrom = (
            f"chr{row['chr']}"
            if not str(row["chr"]).startswith("chr")
            else str(row["chr"])
        )
        pos, ref, alt = int(row["pos"]), str(row["ref"]), str(row["alt"])
        start, end = pos - half_window, pos + half_window
        if start < 1:
            skipped += 1
            continue
        try:
            seq = genome_hg38[chrom][start:end].seq.upper()
        except Exception:
            skipped += 1
            continue
        if len(seq) < 2 * half_window:
            skipped += 1
            continue
        center = half_window
        seq_ref = seq[:center] + ref + seq[center + 1 :]
        seq_alt = seq[:center] + alt + seq[center + 1 :]
        rows.append(
            {
                "variant": row["variant"],
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "pval": row["pval"],
                "sequence_ref": seq_ref,
                "sequence_alt": seq_alt,
            }
        )

    seq_df = pd.DataFrame(rows)
    seq_df.to_csv(out_path, sep="\t", index=False)
    print(f"{len(seq_df):,} sequences extracted (skipped {skipped}) -> {out_path}")
    return seq_df
