from pathlib import Path

import numpy as np
import pandas as pd
from load_counts import (
    DATA_DIR,
    DAY1_ORIGIN,
    N_REPLICATES,
    POPULATIONS,
    RESULTS_DIR,
    SEXES,
    ValidationError,
)

GROUP_KEYS = ["population", "replicate", "sex"]
TOP_FRACTION = 0.10


def top_fraction_mean(ages: pd.Series, fraction: float = TOP_FRACTION) -> float:

    n_top = max(1, int(np.ceil(len(ages) * fraction)))
    return float(np.sort(ages.to_numpy())[-n_top:].mean())


def summarise_replicates(flies: pd.DataFrame) -> pd.DataFrame:

    grouped = flies.groupby(GROUP_KEYS, sort=True)["age_at_death"]

    summary = grouped.agg(
        n="size",
        mean="mean",
        median="median",
        sd="std",
        youngest="min",
        oldest="max",
    )
    summary["top10_mean"] = grouped.apply(top_fraction_mean)

    return summary.reset_index().round(2)


def summarise_populations(flies: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:

    grouped = flies.groupby(["population", "sex"], sort=True)["age_at_death"]

    populations = grouped.agg(
        n_flies="size",
        mean="mean",
        sd="std",
        median="median",
        youngest="min",
        oldest="max",
    )
    # Mean of the longest-lived TOP_FRACTION, the usual stand-in for maximum
    # lifespan: the single oldest fly is one draw and moves with luck.
    populations["top10_mean"] = grouped.apply(top_fraction_mean)

    by_replicate = summary.groupby(["population", "sex"], sort=True)["mean"]
    populations["n_replicates"] = by_replicate.size()
    populations["mean_of_replicate_means"] = by_replicate.mean()
    populations["sd_between_replicates"] = by_replicate.std()

    return populations.reset_index().round(2)


def summarise_spread(summary: pd.DataFrame) -> pd.DataFrame:

    grouped = summary.groupby(["population", "sex"], sort=True)["mean"]

    spread = grouped.agg(
        n_replicates="size",
        mean_of_means="mean",
        sd_between="std",
        lowest="min",
        highest="max",
    )
    spread["range"] = spread["highest"] - spread["lowest"]
    spread["cv"] = spread["sd_between"] / spread["mean_of_means"]

    return spread.reset_index().round(3)


def validate_summary(flies: pd.DataFrame, summary: pd.DataFrame) -> None:
    errors = []

    expected_rows = len(POPULATIONS) * N_REPLICATES * len(SEXES)
    if len(summary) != expected_rows:
        errors.append(f"expected {expected_rows} summary rows, found {len(summary)}")

    if summary["n"].sum() != len(flies):
        errors.append(
            f"summary covers {summary['n'].sum()} flies, the input has {len(flies)}"
        )

    if not (summary["youngest"] <= summary["mean"]).all():
        errors.append("some replicates have a mean below their youngest fly")
    if not (summary["mean"] <= summary["oldest"]).all():
        errors.append("some replicates have a mean above their oldest fly")
    if not (summary["top10_mean"] >= summary["mean"]).all():
        errors.append("some replicates have a top-10% mean below the overall mean")
    if not (summary["top10_mean"] <= summary["oldest"]).all():
        errors.append("some replicates have a top-10% mean above their oldest fly")

    if summary["sd"].isna().any():
        errors.append("some replicates have an undefined standard deviation")

    if errors:
        raise ValidationError(
            "summary failed validation:\n  - " + "\n  - ".join(errors)
        )


def main(flies_csv: Path, out_dir: Path) -> None:
    flies = pd.read_csv(flies_csv)

    summary = summarise_replicates(flies)
    validate_summary(flies, summary)
    populations = summarise_populations(flies, summary)
    spread = summarise_spread(summary)

    out_dir.mkdir(parents=True, exist_ok=True)
    populations.to_csv(out_dir / "lifespan_summary.csv", index=False)
    summary.to_csv(out_dir / "replicate_summary.csv", index=False)
    spread.to_csv(out_dir / "replicate_spread.csv", index=False)

    print(f"day 1 origin: {DAY1_ORIGIN} (assumed, unconfirmed)")
    print("\nlifespan by population (days):")
    print(
        populations.set_index(["population", "sex"])[
            ["n_flies", "mean", "sd", "median", "youngest", "oldest", "top10_mean"]
        ]
        .reindex(POPULATIONS, level="population")
        .to_string()
    )
    print("\nreplicate means (days):")
    print(
        summary.pivot_table(
            index=["population", "sex"], columns="replicate", values="mean"
        )
        .reindex(POPULATIONS, level="population")
        .to_string()
    )
    print("\nbetween-replicate spread:")
    print(
        spread.set_index(["population", "sex"])[
            ["mean_of_means", "sd_between", "range", "cv"]
        ]
        .reindex(POPULATIONS, level="population")
        .to_string()
    )
    print(f"\nwrote 3 tables to {out_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "flies_csv",
        type=Path,
        nargs="?",
        default=DATA_DIR / "flies.csv",
        help="one-row-per-fly table from stage 2",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=RESULTS_DIR,
        help="directory for the summary tables",
    )
    args = parser.parse_args()
    main(args.flies_csv, args.out_dir)
