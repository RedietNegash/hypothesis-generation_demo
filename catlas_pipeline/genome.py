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
