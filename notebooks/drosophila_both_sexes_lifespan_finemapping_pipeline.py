import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
STUDY_URL = "https://dgrpool.epfl.ch/studies/1/get_file?name=summary.tsv"
DEFAULT_OUTPUT_DIR = REPO_DIR / "data" / "finemap" / "both_sexes"
PHENOTYPE_NAME = "longevity_both_sexes_days"
GENOTYPE_URL_PREFIX = "https://zenodo.org/records/837947/files/dgrp2_dm6_dbSNP.vcf"
DEFAULT_GENOTYPE_BFILE = REPO_DIR / "data" / "reference" / "dgrp2_dm6_dbSNP.vcf"
GENOTYPE_MD5 = {
    ".bed": "6e3b0c5b2c186c3a6dec99c20c24e3f4",
    ".bim": "18ba9e60111dd6398d160f9bba89023c",
    ".fam": "d5cbeec02da5e40a89f3c3a9ab8f15cf",
}
FLY_CHROMS = ("2L", "2R", "3L", "3R", "4", "X")
SOURCE_CHROMS = {"2L": "chr2L", "2R": "chr2R", "3L": "chr3L", "3R": "chr3R", "4": "4", "X": "23"}


def download_phenotype_source(destination: Path, url: str = STUDY_URL) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix="arya_study1_", suffix=".tmp", dir=destination.parent, delete=False
    ) as temporary:
        temporary_path = Path(temporary.name)
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                shutil.copyfileobj(response, temporary)
        except BaseException:
            temporary_path.unlink(missing_ok=True)
            raise
    try:
        if temporary_path.stat().st_size == 0:
            raise ValueError(f"Downloaded Arya phenotype source is empty: {url}")
        temporary_path.replace(destination)
    finally:
        temporary_path.unlink(missing_ok=True)
    return destination


def prepare_phenotype(
    output_dir: Path, source_file: Path | None = None,
    source_url: str = STUDY_URL, force_download: bool = False
) -> Path:
    if source_file is None:
        source_file = output_dir / "raw" / "arya_study1_summary.tsv"
        if force_download or not source_file.is_file():
            download_phenotype_source(source_file, source_url)
    elif not source_file.is_file():
        raise FileNotFoundError(f"Arya phenotype source not found: {source_file}")

    source = pd.read_csv(source_file, sep="\t", dtype={"DGRP": "string", "sex": "string"})
    required = {"DGRP", "sex", "mn_Longevity"}
    if not required.issubset(source.columns):
        raise ValueError(f"Arya phenotype source is missing columns: {sorted(required - set(source.columns))}")
    source = source.loc[source["sex"].isin(["F", "M"]), list(required)].copy()
    if source.empty or source[["DGRP", "sex"]].isna().any().any():
        raise ValueError("Arya phenotype source has no valid male and female line identifiers")
    if source.duplicated(["DGRP", "sex"]).any():
        raise ValueError("Arya phenotype source contains duplicate line-sex measurements")
    source["mn_Longevity"] = pd.to_numeric(source["mn_Longevity"], errors="coerce")
    source.loc[~np.isfinite(source["mn_Longevity"]), "mn_Longevity"] = np.nan
    paired = source.pivot(index="DGRP", columns="sex", values="mn_Longevity")
    if not {"F", "M"}.issubset(paired.columns):
        raise ValueError("Arya phenotype source must contain both male and female longevity")
    paired = paired.dropna(subset=["F", "M"]).copy()
    if paired.empty or (paired[["F", "M"]] <= 0).any().any():
        raise ValueError("No paired, positive male and female longevity values found")

    phenotype = paired.rename(columns={"F": "female_days", "M": "male_days"}).reset_index()
    identifiers = phenotype["DGRP"].str.extract(r"^DGRP_(\d+)$", expand=False)
    if identifiers.isna().any():
        raise ValueError("Arya phenotype source contains malformed DGRP line identifiers")
    phenotype.insert(1, "FID", "line")
    phenotype.insert(2, "IID", identifiers.astype(int).astype(str))
    if phenotype["IID"].duplicated().any():
        raise ValueError("Arya phenotype source contains duplicate normalized line identifiers")
    phenotype[PHENOTYPE_NAME] = phenotype[["female_days", "male_days"]].mean(axis=1)
    phenotype = phenotype[
        ["DGRP", "FID", "IID", "female_days", "male_days", PHENOTYPE_NAME]
    ].sort_values("IID", key=lambda values: values.astype(int))

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "arya_both_sexes_phenotype.tsv"
    temporary_file = output_file.with_suffix(".tsv.tmp")
    try:
        phenotype.to_csv(temporary_file, sep="\t", index=False)
        temporary_file.replace(output_file)
    finally:
        temporary_file.unlink(missing_ok=True)
    digest = hashlib.sha256(source_file.read_bytes()).hexdigest()
    print(
        f"Prepared {len(phenotype)} paired DGRP lines from {source_file}; "
        f"source SHA-256 {digest}: {output_file}"
    )
    return output_file


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_genotype_reference(prefix: Path, download: bool) -> None:
    for extension, expected_md5 in GENOTYPE_MD5.items():
        destination = Path(f"{prefix}{extension}")
        if not download:
            if not destination.is_file():
                raise FileNotFoundError(f"DGRP2 genotype reference is missing: {destination}")
            continue
        if destination.is_file() and file_digest(destination, "md5") == expected_md5:
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            prefix="dgrp2_", suffix=".tmp", dir=destination.parent, delete=False
        ) as temporary:
            temporary_path = Path(temporary.name)
            try:
                with urllib.request.urlopen(
                    GENOTYPE_URL_PREFIX + extension + "?download=1", timeout=120
                ) as response:
                    shutil.copyfileobj(response, temporary)
            except BaseException:
                temporary_path.unlink(missing_ok=True)
                raise
        try:
            observed_md5 = file_digest(temporary_path, "md5")
            if observed_md5 != expected_md5:
                raise ValueError(
                    f"DGRP2 genotype checksum mismatch for {extension}: {observed_md5}"
                )
            temporary_path.replace(destination)
        finally:
            temporary_path.unlink(missing_ok=True)


def prepare_genotype_reference(
    phenotype_file: Path, bfile_prefix: Path, output_dir: Path,
    download: bool = True
) -> Path:
    if not phenotype_file.is_file():
        raise FileNotFoundError(f"Prepared Arya phenotype not found: {phenotype_file}")
    ensure_genotype_reference(bfile_prefix, download)
    phenotype = pd.read_csv(phenotype_file, sep="\t", dtype={"FID": "string", "IID": "string"})
    if not {"FID", "IID", PHENOTYPE_NAME}.issubset(phenotype.columns):
        raise ValueError("Prepared Arya phenotype lacks FID, IID, or mean longevity")
    if phenotype[["FID", "IID"]].isna().any().any() or phenotype.duplicated(
        ["FID", "IID"]
    ).any():
        raise ValueError("Prepared Arya phenotype has missing or duplicate sample IDs")
    fam = pd.read_csv(
        Path(f"{bfile_prefix}.fam"), sep=r"\s+", header=None,
        usecols=[0, 1], names=["FID", "IID"], dtype="string",
    )
    if fam.empty or fam.isna().any().any() or fam.duplicated(["FID", "IID"]).any():
        raise ValueError("DGRP2 reference has missing or duplicate sample IDs")
    matched = phenotype[["FID", "IID"]].merge(
        fam, on=["FID", "IID"], how="inner", validate="one_to_one"
    )
    if len(matched) < 2:
        raise ValueError("Fewer than two paired Arya lines match the DGRP2 genotypes")
    output_dir.mkdir(parents=True, exist_ok=True)
    keep_file = output_dir / "arya_genotyped_lines.keep"
    temporary_file = keep_file.with_suffix(".keep.tmp")
    try:
        matched.to_csv(temporary_file, sep="\t", header=False, index=False)
        temporary_file.replace(keep_file)
    finally:
        temporary_file.unlink(missing_ok=True)
    print(
        f"Matched {len(matched)} of {len(phenotype)} paired Arya lines "
        f"to DGRP2 genotypes: {keep_file}"
    )
    return keep_file


def plink_files_exist(prefix: Path) -> bool:
    return all(Path(f"{prefix}{extension}").is_file() for extension in (".bed", ".bim", ".fam"))


def run_command(arguments: list[str]) -> None:
    print("Running:", " ".join(arguments), flush=True)
    subprocess.run(arguments, check=True)


def prepare_qc_and_pca(
    bfile_prefix: Path, keep_file: Path, output_dir: Path,
    plink_bin: str = "plink", plink2_bin: str = "plink2", force: bool = False,
) -> Path:
    if not keep_file.is_file() or keep_file.stat().st_size == 0:
        raise FileNotFoundError(f"Matched Arya lines are missing: {keep_file}")
    if not plink_files_exist(bfile_prefix):
        raise FileNotFoundError(f"DGRP2 genotype reference is incomplete: {bfile_prefix}")
    plink = shutil.which(plink_bin)
    plink2 = shutil.which(plink2_bin)
    if plink is None or plink2 is None:
        raise FileNotFoundError("PLINK 1.9 and PLINK 2 must be installed or passed by path")

    qc_dir = output_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    for chrom in FLY_CHROMS:
        prefix = qc_dir / chrom
        if plink_files_exist(prefix) and not force:
            continue
        run_command([
            plink2, "--bfile", str(bfile_prefix), "--chr", SOURCE_CHROMS[chrom],
            "--keep", str(keep_file), "--maf", "0.01", "--geno", "0.05",
            "--allow-extra-chr", "--make-bed", "--out", str(prefix),
        ])
        if not plink_files_exist(prefix):
            raise RuntimeError(f"PLINK did not create chromosome {chrom}: {prefix}")
        bim = Path(f"{prefix}.bim")
        temporary_bim = Path(f"{bim}.tmp")
        try:
            with bim.open() as source, temporary_bim.open("w") as destination:
                for line in source:
                    fields = line.rstrip("\n").split("\t", 1)
                    if len(fields) != 2 or fields[0] not in (SOURCE_CHROMS[chrom], chrom):
                        raise ValueError(f"Unexpected chromosome label in {bim}: {line.strip()}")
                    destination.write(f"{chrom}\t{fields[1]}\n")
            temporary_bim.replace(bim)
        finally:
            temporary_bim.unlink(missing_ok=True)

    merged = qc_dir / "merged_qc"
    if not plink_files_exist(merged) or force:
        merge_list = qc_dir / "merge_list.txt"
        merge_list.write_text("".join(
            f"{qc_dir / chrom}.bed {qc_dir / chrom}.bim {qc_dir / chrom}.fam\n"
            for chrom in FLY_CHROMS[1:]
        ))
        run_command([
            plink, "--bfile", str(qc_dir / FLY_CHROMS[0]), "--merge-list",
            str(merge_list), "--allow-extra-chr", "--make-bed", "--out", str(merged),
        ])
        if not plink_files_exist(merged):
            raise RuntimeError(f"PLINK did not create the merged QC reference: {merged}")

    pca_dir = output_dir / "pca"
    pca_dir.mkdir(parents=True, exist_ok=True)
    prune_prefix = pca_dir / "ld_pruned"
    prune_file = Path(f"{prune_prefix}.prune.in")
    if not prune_file.is_file() or force:
        run_command([
            plink2, "--bfile", str(merged), "--allow-extra-chr", "--indep-pairwise",
            "200kb", "1", "0.2", "--out", str(prune_prefix),
        ])
        if not prune_file.is_file() or prune_file.stat().st_size == 0:
            raise RuntimeError("PLINK produced no variants for PCA after LD pruning")
    pca_prefix = pca_dir / "arya_pca"
    eigenvec = Path(f"{pca_prefix}.eigenvec")
    if not eigenvec.is_file() or force:
        run_command([
            plink2, "--bfile", str(merged), "--extract", str(prune_file),
            "--allow-extra-chr", "--pca", "10", "--out", str(pca_prefix),
        ])
        if not eigenvec.is_file():
            raise RuntimeError(f"PLINK did not create PCA covariates: {eigenvec}")
    print(f"Prepared QC genotypes and 10 PCs: {eigenvec}")
    return eigenvec


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Arya et al. DGRP male and female longevity fine-mapping pipeline"
    )
    parser.add_argument("--stage", choices=("phenotype", "genotype", "qc"), default="phenotype")
    parser.add_argument("--phenotype-source", type=Path)
    parser.add_argument("--phenotype-url", default=STUDY_URL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--genotype-bfile", type=Path)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plink-bin", default="plink")
    parser.add_argument("--plink2-bin", default="plink2")
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser()
    if args.stage == "phenotype":
        prepare_phenotype(
            output_dir,
            args.phenotype_source.expanduser() if args.phenotype_source else None,
            args.phenotype_url, args.force_download,
        )
    elif args.stage == "genotype":
        prepare_genotype_reference(
            output_dir / "arya_both_sexes_phenotype.tsv",
            args.genotype_bfile.expanduser() if args.genotype_bfile else DEFAULT_GENOTYPE_BFILE,
            output_dir, download=args.genotype_bfile is None,
        )
    elif args.stage == "qc":
        prepare_qc_and_pca(
            args.genotype_bfile.expanduser() if args.genotype_bfile else DEFAULT_GENOTYPE_BFILE,
            output_dir / "arya_genotyped_lines.keep", output_dir,
            plink_bin=args.plink_bin, plink2_bin=args.plink2_bin, force=args.force,
        )


if __name__ == "__main__":
    main()
