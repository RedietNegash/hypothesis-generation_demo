from pathlib import Path

import numpy as np
import pandas as pd
from load_counts import (
    AGE_BASE,
    AGE_OFFSET_DAYS,
    DATA_DIR,
    DAY1_ORIGIN,
    EXPECTED_N,
    N_REPLICATES,
    POPULATIONS,
    SEXES,
    ValidationError,
)

GROUP_KEYS = ["population", "replicate", "sex"]


def day_to_age(day: "pd.Series | np.ndarray") -> np.ndarray:

    return np.asarray(day) - AGE_BASE + 1 + AGE_OFFSET_DAYS


def expand(deaths: pd.DataFrame) -> pd.DataFrame:

    counts = deaths["deaths"].to_numpy()

    flies = pd.DataFrame(
        {
            "population": np.repeat(deaths["population"].to_numpy(), counts),
            "replicate": np.repeat(deaths["replicate"].to_numpy(), counts),
            "sex": np.repeat(deaths["sex"].to_numpy(), counts),
            "age_at_death": np.repeat(day_to_age(deaths["day"]), counts),
        }
    )
    flies["event_observed"] = 1

    return flies


def validate_expansion(deaths: pd.DataFrame, flies: pd.DataFrame) -> None:

    errors = []

    expected_rows = int(deaths["deaths"].sum())
    if len(flies) != expected_rows:
        errors.append(f"expanded to {len(flies)} rows, counts sum to {expected_rows}")

    sizes = flies.groupby(GROUP_KEYS, sort=True).size()
    bad_sizes = sizes[sizes != EXPECTED_N]
    for (pop, rep, sex), n in bad_sizes.items():
        errors.append(f"{pop} replicate {rep} {sex}: {n} flies, expected {EXPECTED_N}")

    expected_groups = len(POPULATIONS) * N_REPLICATES * len(SEXES)
    if len(sizes) != expected_groups:
        errors.append(
            f"expected {expected_groups} groups after expansion, found {len(sizes)}"
        )

    rebuilt = (
        flies.groupby(GROUP_KEYS + ["age_at_death"], sort=True)
        .size()
        .rename("deaths")
        .reset_index()
    )
    source = deaths[deaths["deaths"] > 0].copy()
    source["age_at_death"] = day_to_age(source["day"])
    source = source[GROUP_KEYS + ["age_at_death", "deaths"]].sort_values(
        GROUP_KEYS + ["age_at_death"]
    )

    merged = source.merge(
        rebuilt,
        on=GROUP_KEYS + ["age_at_death"],
        how="outer",
        suffixes=("_source", "_rebuilt"),
        indicator=True,
    )
    mismatched = merged[
        (merged["_merge"] != "both")
        | (merged["deaths_source"] != merged["deaths_rebuilt"])
    ]
    if not mismatched.empty:
        errors.append(
            f"{len(mismatched)} age/count cells do not round-trip; "
            f"first: {mismatched.iloc[0].to_dict()}"
        )

    if errors:
        raise ValidationError(
            "expansion failed validation:\n  - " + "\n  - ".join(errors)
        )


def summarise(flies: pd.DataFrame) -> pd.DataFrame:
    return (
        flies.groupby("population", sort=False)["age_at_death"]
        .agg(flies="size", mean="mean", median="median", oldest="max")
        .round(1)
        .reindex(POPULATIONS)
    )


def main(deaths_csv: Path, output: Path) -> None:
    deaths = pd.read_csv(deaths_csv)
    flies = expand(deaths)
    validate_expansion(deaths, flies)

    output.parent.mkdir(parents=True, exist_ok=True)
    flies.to_csv(output, index=False)

    print(f"day 1 origin: {DAY1_ORIGIN} (assumed, unconfirmed)")
    print(f"age_at_death = day - {AGE_BASE} + 1 + {AGE_OFFSET_DAYS}")
    print(summarise(flies).to_string())
    print(f"\nwrote {len(flies)} rows to {output}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "deaths_csv",
        type=Path,
        nargs="?",
        default=DATA_DIR / "deaths_long.csv",
        help="long-form death counts from stage 1",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DATA_DIR / "flies.csv",
        help="where to write the one-row-per-fly table",
    )
    args = parser.parse_args()
    main(args.deaths_csv, args.output)
