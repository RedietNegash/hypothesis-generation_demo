from pathlib import Path

import pandas as pd

# --- Configuration ----------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"

DAY1_ORIGIN = "eclosion"

AGE_BASE = 1


AGE_OFFSET_DAYS = 0

POPULATIONS = ("B", "O", "SO")
SHEET_TEMPLATE = "{population} death"
N_REPLICATES = 5
SEXES = ("Males", "Females")
EXPECTED_N = 250  # flies per replicate per sex at the start of the assay


class ValidationError(Exception):
    pass


# --- Loading ----------------------------------------------------------------


def load_death_sheet(workbook: Path, population: str) -> pd.DataFrame:

    sheet = SHEET_TEMPLATE.format(population=population)
    wide = pd.read_excel(workbook, sheet_name=sheet)

    expected_cols = 2 + N_REPLICATES * len(SEXES)
    if wide.shape[1] != expected_cols:
        raise ValidationError(
            f"{sheet}: expected {expected_cols} columns "
            f"(Date, Day, then {N_REPLICATES} Males/Females pairs), "
            f"found {wide.shape[1]}"
        )

    count_cols = wide.columns[2:]
    frames = []
    for offset, col in enumerate(count_cols):
        replicate = offset // len(SEXES) + 1
        sex = SEXES[offset % len(SEXES)]

        if not str(col).startswith(sex):
            raise ValidationError(
                f"{sheet}: column {offset + 3} is '{col}', expected a "
                f"'{sex}' column for replicate {replicate}. The Males/Females "
                f"pairs are not in the order this loader assumes."
            )
        counts = wide[col]
        last_recorded = counts.last_valid_index()
        if last_recorded is not None:
            interior_blanks = counts.iloc[: last_recorded + 1].isna()
            if interior_blanks.any():
                days = wide.loc[interior_blanks[interior_blanks].index, "Day"]
                raise ValidationError(
                    f"{sheet}: column '{col}' is blank on day(s) "
                    f"{list(days)} but has later entries. A blank between "
                    f"recorded deaths is a missed census, not zero deaths."
                )

        frames.append(
            pd.DataFrame(
                {
                    "population": population,
                    "replicate": replicate,
                    "sex": sex,
                    "day": wide["Day"].astype(int),
                    "deaths": counts.fillna(0).astype(int),
                }
            )
        )

    return pd.concat(frames, ignore_index=True)


def load_all(workbook: Path) -> pd.DataFrame:
    return pd.concat(
        [load_death_sheet(workbook, pop) for pop in POPULATIONS],
        ignore_index=True,
    )


# --- Validation -------------------------------------------------------------


def validate(deaths: pd.DataFrame) -> None:

    errors = []

    if (deaths["deaths"] < 0).any():
        errors.append("negative death counts present")

    groups = deaths.groupby(["population", "replicate", "sex"], sort=True)

    expected_groups = len(POPULATIONS) * N_REPLICATES * len(SEXES)
    if groups.ngroups != expected_groups:
        errors.append(
            f"expected {expected_groups} population/replicate/sex groups, "
            f"found {groups.ngroups}"
        )

    totals = groups["deaths"].sum()
    bad_totals = totals[totals != EXPECTED_N]
    for (pop, rep, sex), total in bad_totals.items():
        errors.append(
            f"{pop} replicate {rep} {sex}: deaths sum to {total}, "
            f"expected {EXPECTED_N} (difference {total - EXPECTED_N:+d})"
        )

    for (pop, rep, sex), group in groups:
        days = group["day"].to_numpy()
        expected_days = range(1, len(days) + 1)
        if list(days) != list(expected_days):
            errors.append(
                f"{pop} replicate {rep} {sex}: days are not consecutive "
                f"from 1 (first={days[0]}, last={days[-1]}, n={len(days)})"
            )

    if errors:
        raise ValidationError(
            "workbook failed validation:\n  - " + "\n  - ".join(errors)
        )


def summarise(deaths: pd.DataFrame) -> pd.DataFrame:
    return (
        deaths.groupby("population", sort=False)
        .agg(
            replicates=("replicate", "nunique"),
            last_day=("day", "max"),
            total_deaths=("deaths", "sum"),
        )
        .reindex(POPULATIONS)
    )


# --- Entry point ------------------------------------------------------------


def main(workbook: Path, output: Path) -> None:
    deaths = load_all(workbook)
    validate(deaths)

    output.parent.mkdir(parents=True, exist_ok=True)
    deaths.to_csv(output, index=False)

    print(f"day 1 origin: {DAY1_ORIGIN} (assumed, unconfirmed)")
    print(summarise(deaths).to_string())
    print(f"\nwrote {len(deaths)} rows to {output}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path, help="path to the counts .xlsx")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DATA_DIR / "deaths_long.csv",
        help="where to write the long-form death counts",
    )
    args = parser.parse_args()
    main(args.workbook, args.output)
