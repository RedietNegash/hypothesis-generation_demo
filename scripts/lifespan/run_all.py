import subprocess
import sys
from pathlib import Path

from load_counts import BASE_DIR, DATA_DIR, RESULTS_DIR


STAGES = [
    ("load_counts.py", "read the workbook into long-form death counts"),
    ("expand_ages.py", "expand counts into one row per fly"),
    ("summarise_replicates.py", "lifespan summaries by population and replicate"),
    ("survival_curves.py", "Kaplan-Meier survival curves and death rates"),
    ("compare_populations.py", "compare B, O and SO across replicates"),
    ("plot_survival.py", "render the survival curve figures"),
]

DEFAULT_WORKBOOK = BASE_DIR / "counts-8nqxycmkc7ysxro9d1ucfisuio.xlsx"


def run_stage(script: str, arguments: list[str]) -> None:
    command = [sys.executable, str(BASE_DIR / script), *arguments]
    result = subprocess.run(command, cwd=BASE_DIR)
    if result.returncode != 0:
        raise SystemExit(
            f"\nstage '{script}' exited with code {result.returncode}; "
            f"later stages were not run"
        )


def main(workbook: Path) -> None:
    if not workbook.is_file():
        raise SystemExit(f"workbook not found: {workbook}")

    for number, (script, description) in enumerate(STAGES, start=1):
        banner = f" stage {number}/{len(STAGES)}: {description} "
        print(f"\n{banner:=^78}\n")
        run_stage(script, [str(workbook)] if number == 1 else [])

    print(f"\n{' done ':=^78}\n")
    print(f"intermediate tables: {DATA_DIR}")
    print(f"results:             {RESULTS_DIR}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "workbook",
        type=Path,
        nargs="?",
        default=DEFAULT_WORKBOOK,
        help="path to the counts .xlsx",
    )
    args = parser.parse_args()
    main(args.workbook)
