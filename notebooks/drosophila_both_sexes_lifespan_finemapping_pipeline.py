import argparse
import hashlib
from pathlib import Path
import shutil
import tempfile
import urllib.request

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
STUDY_URL = "https://dgrpool.epfl.ch/studies/1/get_file?name=summary.tsv"
DEFAULT_OUTPUT_DIR = REPO_DIR / "data" / "finemap" / "both_sexes"
PHENOTYPE_NAME = "longevity_both_sexes_days"


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Arya et al. DGRP male and female longevity fine-mapping pipeline"
    )
    parser.add_argument("--stage", choices=("phenotype",), default="phenotype")
    parser.add_argument("--phenotype-source", type=Path)
    parser.add_argument("--phenotype-url", default=STUDY_URL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()
    prepare_phenotype(
        args.output_dir.expanduser(),
        args.phenotype_source.expanduser() if args.phenotype_source else None,
        args.phenotype_url, args.force_download,
    )


if __name__ == "__main__":
    main()
