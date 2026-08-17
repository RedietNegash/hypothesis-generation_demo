from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from load_counts import (
    AGE_OFFSET_DAYS,
    DAY1_ORIGIN,
    DATA_DIR,
    POPULATIONS,
    RESULTS_DIR,
    ValidationError,
)

GROUP_KEYS = ["population", "replicate", "sex"]
POOLED_KEYS = ["population", "sex"]


MORTALITY_BIN_DAYS = 7


MIN_AT_RISK = 20


def _km_frame(ages: pd.Series, events: pd.Series, label: str) -> pd.DataFrame:
    kmf = KaplanMeierFitter()
    kmf.fit(ages, event_observed=events, label=label)

    curve = kmf.survival_function_.copy()
    curve.columns = ["survival"]
    ci = kmf.confidence_interval_
    curve["ci_lower"] = ci.iloc[:, 0].to_numpy()
    curve["ci_upper"] = ci.iloc[:, 1].to_numpy()
    curve["at_risk"] = [int((ages >= t).sum()) for t in curve.index]

    curve = curve.reset_index().rename(columns={"timeline": "age"})
    curve["median_lifespan"] = kmf.median_survival_time_
    return curve


def survival_by(flies: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    frames = []
    for values, group in flies.groupby(keys, sort=True):
        values = values if isinstance(values, tuple) else (values,)
        curve = _km_frame(
            group["age_at_death"],
            group["event_observed"],
            label="/".join(str(v) for v in values),
        )
        for key, value in zip(keys, values):
            curve.insert(0, key, value)
        frames.append(curve.drop(columns="label", errors="ignore"))

    return pd.concat(frames, ignore_index=True)


def mortality_by(
    flies: pd.DataFrame,
    keys: list[str],
    bin_days: int = MORTALITY_BIN_DAYS,
    min_at_risk: int = MIN_AT_RISK,
) -> pd.DataFrame:

    origin = AGE_OFFSET_DAYS  # cohort age when the assay began

    rows = []
    for values, group in flies.groupby(keys, sort=True):
        values = values if isinstance(values, tuple) else (values,)
        ages = group["age_at_death"].to_numpy()
        n_total = len(ages)

        edges = np.arange(origin, ages.max() + bin_days, bin_days)
        for start in edges:
            end = start + bin_days
            at_risk = int((ages > start).sum())
            if at_risk < min_at_risk:
                continue
            deaths = int(((ages > start) & (ages <= end)).sum())
            q = deaths / at_risk

            # q == 1 wipes out the whole risk set, and -ln(0) is infinite:
            # the data fix no upper bound on the rate.
            hazard = -np.log(1 - q) / bin_days if q < 1 else np.nan

            rows.append(
                dict(
                    zip(keys, values),
                    bin_start=int(start),
                    bin_end=int(end),
                    bin_midpoint=start + bin_days / 2,
                    at_risk=at_risk,
                    deaths=deaths,
                    death_probability=q,
                    hazard_rate=hazard,
                    cohort_n=n_total,
                )
            )

    return pd.DataFrame(rows).round(6)


def validate_survival(flies: pd.DataFrame, pooled: pd.DataFrame) -> None:
    errors = []

    if not pooled["survival"].between(0, 1).all():
        errors.append("survival values outside [0, 1]")

    for values, curve in pooled.groupby(POOLED_KEYS, sort=True):
        label = "/".join(str(v) for v in values)
        s = curve.sort_values("age")["survival"].to_numpy()

        if not np.all(np.diff(s) <= 1e-12):
            errors.append(f"{label}: survival is not monotonically decreasing")
        if s[-1] > 1e-9:
            errors.append(
                f"{label}: survival ends at {s[-1]:.4f}, expected 0 "
                f"(no censoring means every fly should be dead)"
            )

    for values, group in flies.groupby(POOLED_KEYS, sort=True):
        label = "/".join(str(v) for v in values)
        km_median = pooled.loc[
            (pooled["population"] == values[0]) & (pooled["sex"] == values[1]),
            "median_lifespan",
        ].iloc[0]
        raw_median = group["age_at_death"].median()
        if abs(km_median - raw_median) > 1.0:
            errors.append(f"{label}: KM median {km_median} vs empirical {raw_median}")

    if errors:
        raise ValidationError(
            "survival curves failed validation:\n  - " + "\n  - ".join(errors)
        )


def main(flies_csv: Path, out_dir: Path) -> None:
    flies = pd.read_csv(flies_csv)

    per_replicate = survival_by(flies, GROUP_KEYS)
    pooled = survival_by(flies, POOLED_KEYS)
    validate_survival(flies, pooled)

    mortality = mortality_by(flies, POOLED_KEYS)
    mortality_replicate = mortality_by(flies, GROUP_KEYS)

    out_dir.mkdir(parents=True, exist_ok=True)
    per_replicate.to_csv(out_dir / "survival_curves.csv", index=False)
    pooled.to_csv(out_dir / "survival_pooled.csv", index=False)
    mortality.to_csv(out_dir / "mortality.csv", index=False)
    mortality_replicate.to_csv(out_dir / "mortality_replicate.csv", index=False)

    print(f"day 1 origin: {DAY1_ORIGIN} (assumed, unconfirmed)")

    medians = (
        pooled.groupby(POOLED_KEYS)["median_lifespan"]
        .first()
        .unstack()
        .reindex(POPULATIONS)
    )
    print("\nKaplan-Meier median lifespan (days):")
    print(medians.to_string())

    print("\nage at which survival first drops below 10%:")
    tail = (
        pooled[pooled["survival"] < 0.10]
        .groupby(POOLED_KEYS)["age"]
        .min()
        .unstack()
        .reindex(POPULATIONS)
    )
    print(tail.to_string())

    print(
        f"\nearly deaths, first two {MORTALITY_BIN_DAYS}-day bins "
        f"(fraction of flies alive at bin start):"
    )
    early = (
        mortality[mortality["bin_start"] < 2 * MORTALITY_BIN_DAYS]
        .pivot_table(index=POOLED_KEYS, columns="bin_start", values="death_probability")
        .reindex(POPULATIONS, level="population")
        .round(4)
    )
    print(early.to_string())

    print(f"\nwrote 4 tables to {out_dir}")


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
        help="directory for the survival and mortality tables",
    )
    args = parser.parse_args()
    main(args.flies_csv, args.out_dir)
