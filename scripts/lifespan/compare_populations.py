from pathlib import Path

import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test
from load_counts import (
    DATA_DIR,
    N_REPLICATES,
    POPULATIONS,
    RESULTS_DIR,
    ValidationError,
)
from scipy import stats

GROUP_KEYS = ["population", "replicate", "sex"]
POPULATION_PAIRS = [("B", "O"), ("O", "SO"), ("B", "SO")]
MEASURES = ["mean_lifespan", "median_lifespan", "top10_mean"]
TOP_FRACTION = 0.10

# p-values reach 1e-14; rounding to fixed decimals would collapse them to zero.
UNROUNDED = ["p_value"]


def round_except_p(frame: pd.DataFrame, decimals: int = 6) -> pd.DataFrame:
    return frame.assign(
        **{
            column: frame[column].round(decimals)
            for column in frame.columns
            if column not in UNROUNDED and pd.api.types.is_numeric_dtype(frame[column])
        }
    )


def replicate_measures(flies: pd.DataFrame) -> pd.DataFrame:
    ages = flies.groupby(GROUP_KEYS, sort=True)["age_at_death"]

    measures = ages.agg(n="size", mean_lifespan="mean", median_lifespan="median")
    measures["top10_mean"] = ages.apply(
        lambda a: float(
            np.sort(a.to_numpy())[-max(1, int(np.ceil(len(a) * TOP_FRACTION))) :].mean()
        )
    )
    return measures.reset_index()


def compare_pairs(measures: pd.DataFrame) -> pd.DataFrame:

    rows = []
    for sex, by_sex in measures.groupby("sex", sort=True):
        for measure in MEASURES:
            for low, high in POPULATION_PAIRS:
                a = by_sex.loc[by_sex["population"] == low, measure].to_numpy()
                b = by_sex.loc[by_sex["population"] == high, measure].to_numpy()

                test = stats.ttest_ind(a, b, equal_var=False)
                pooled_sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)

                rows.append(
                    {
                        "sex": sex,
                        "measure": measure,
                        "population_low": low,
                        "population_high": high,
                        "n_low": len(a),
                        "n_high": len(b),
                        "mean_low": a.mean(),
                        "mean_high": b.mean(),
                        "difference": b.mean() - a.mean(),
                        "ratio": b.mean() / a.mean(),
                        "cohens_d": (b.mean() - a.mean()) / pooled_sd,
                        "t_statistic": test.statistic,
                        "p_value": test.pvalue,
                    }
                )

    return round_except_p(pd.DataFrame(rows))


def compare_sexes(measures: pd.DataFrame) -> pd.DataFrame:

    rows = []
    for population, group in measures.groupby("population", sort=True):
        wide = group.pivot(
            index="replicate", columns="sex", values="mean_lifespan"
        ).dropna()
        test = stats.ttest_rel(wide["Females"], wide["Males"])
        difference = (wide["Females"] - wide["Males"]).to_numpy()

        rows.append(
            {
                "population": population,
                "n_replicates": len(wide),
                "mean_females": wide["Females"].mean(),
                "mean_males": wide["Males"].mean(),
                "difference": difference.mean(),
                "sd_of_difference": difference.std(ddof=1),
                "t_statistic": test.statistic,
                "p_value": test.pvalue,
            }
        )

    return round_except_p(pd.DataFrame(rows))


def pooled_logrank(flies: pd.DataFrame) -> pd.DataFrame:

    rows = []
    for sex, by_sex in flies.groupby("sex", sort=True):
        for low, high in POPULATION_PAIRS:
            a = by_sex.loc[by_sex["population"] == low, "age_at_death"]
            b = by_sex.loc[by_sex["population"] == high, "age_at_death"]
            result = logrank_test(a, b, np.ones(len(a)), np.ones(len(b)))

            rows.append(
                {
                    "sex": sex,
                    "population_low": low,
                    "population_high": high,
                    "n_flies_low": len(a),
                    "n_flies_high": len(b),
                    "chi_squared": result.test_statistic,
                    "p_value": result.p_value,
                    "note": "pseudoreplicated; flies are not independent",
                }
            )

    return round_except_p(pd.DataFrame(rows))


def validate_comparison(measures: pd.DataFrame, pairs: pd.DataFrame) -> None:
    errors = []

    expected_rows = len(POPULATIONS) * N_REPLICATES * 2
    if len(measures) != expected_rows:
        errors.append(
            f"expected {expected_rows} replicate/sex rows, found {len(measures)}"
        )

    counts = measures.groupby(["population", "sex"]).size()
    if (counts != N_REPLICATES).any():
        errors.append(
            f"not {N_REPLICATES} replicates per population/sex: "
            f"{dict(counts[counts != N_REPLICATES])}"
        )

    if not pairs["p_value"].between(0, 1).all():
        errors.append("p-values outside [0, 1]")

    if not ((pairs["n_low"] == N_REPLICATES) & (pairs["n_high"] == N_REPLICATES)).all():
        errors.append(
            f"a pairwise test used something other than the {N_REPLICATES} "
            f"replicates as its sample"
        )

    ordering = measures.groupby("population")["mean_lifespan"].mean()
    if not (ordering["B"] < ordering["O"] < ordering["SO"]):
        errors.append(
            f"populations do not rank B < O < SO by mean lifespan: {dict(ordering)}"
        )

    if errors:
        raise ValidationError(
            "population comparison failed validation:\n  - " + "\n  - ".join(errors)
        )


def show(frame: pd.DataFrame, decimals: int = 4) -> str:

    formatters = {
        column: (
            "{:.2e}".format if column in UNROUNDED else f"{{:.{decimals}f}}".format
        )
        for column in frame.columns
        if pd.api.types.is_float_dtype(frame[column])
    }
    return frame.to_string(formatters=formatters)


def main(flies_csv: Path, out_dir: Path) -> None:
    flies = pd.read_csv(flies_csv)

    measures = replicate_measures(flies)
    pairs = compare_pairs(measures)
    validate_comparison(measures, pairs)

    sexes = compare_sexes(measures)
    logrank = pooled_logrank(flies)

    out_dir.mkdir(parents=True, exist_ok=True)
    measures.round(4).to_csv(out_dir / "replicate_measures.csv", index=False)
    pairs.to_csv(out_dir / "population_comparison.csv", index=False)
    sexes.to_csv(out_dir / "sex_comparison.csv", index=False)
    logrank.to_csv(out_dir / "pooled_logrank.csv", index=False)

    print(f"unit of analysis: replicate population (n={N_REPLICATES} per group)")

    for measure in MEASURES:
        print(f"\n{measure}, replicates as the unit (Welch's t-test):")
        print(
            show(
                pairs[pairs["measure"] == measure].set_index(
                    ["sex", "population_low", "population_high"]
                )[
                    [
                        "mean_low",
                        "mean_high",
                        "difference",
                        "ratio",
                        "cohens_d",
                        "p_value",
                    ]
                ]
            )
        )

    print("\nfemale minus male mean lifespan, paired within replicate:")
    print(
        show(
            sexes.set_index("population")[
                ["mean_females", "mean_males", "difference", "p_value"]
            ].reindex(POPULATIONS),
            decimals=3,
        )
    )

    print("\nfor contrast, the same pairs pooled over flies (do not quote):")
    print(
        show(
            logrank.set_index(["sex", "population_low", "population_high"])[
                ["n_flies_low", "n_flies_high", "chi_squared", "p_value"]
            ],
            decimals=1,
        )
    )

    print(f"\nwrote 4 tables to {out_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
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
        help="directory for the comparison tables",
    )
    args = parser.parse_args()
    main(args.flies_csv, args.out_dir)
