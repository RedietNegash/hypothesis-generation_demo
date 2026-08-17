from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from load_counts import POPULATIONS, RESULTS_DIR, SEXES, ValidationError


POPULATION_COLOURS = {"B": "#2a78d6", "O": "#eb6834", "SO": "#1baf7a"}

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"

DPI = 200


def apply_chrome(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis="y", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, length=4, width=0.8)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])


def label_curve_end(ax: plt.Axes, curve: pd.DataFrame, text: str) -> None:

    alive = curve[curve["survival"] > 0]
    tail = alive.iloc[-1]
    ax.annotate(
        text,
        xy=(tail["age"], tail["survival"]),
        xytext=(6, 6),
        textcoords="offset points",
        color=INK,
        fontsize=9,
    )


def plot_populations(pooled: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(
        1, len(SEXES), figsize=(11, 4.4), sharey=True, facecolor=SURFACE
    )

    for ax, sex in zip(axes, SEXES):
        apply_chrome(ax)
        for population in POPULATIONS:
            curve = pooled[
                (pooled["population"] == population) & (pooled["sex"] == sex)
            ].sort_values("age")
            colour = POPULATION_COLOURS[population]

            ax.fill_between(
                curve["age"],
                curve["ci_lower"],
                curve["ci_upper"],
                step="post",
                color=colour,
                alpha=0.15,
                linewidth=0,
                zorder=2,
            )
            ax.step(
                curve["age"],
                curve["survival"],
                where="post",
                color=colour,
                linewidth=2,
                label=population,
                zorder=3,
            )
            label_curve_end(ax, curve, population)

        ax.set_title(sex, color=INK, fontsize=11, loc="left", pad=10)
        ax.set_xlabel("Age from eclosion (days)", color=INK_SECONDARY, fontsize=9)
        ax.set_xlim(0, pooled["age"].max() * 1.06)
        ax.set_ylim(0, 1.02)

    axes[0].set_ylabel("Surviving fraction", color=INK_SECONDARY, fontsize=9)
    axes[0].legend(
        title="Population",
        frameon=False,
        labelcolor=INK_SECONDARY,
        fontsize=9,
        title_fontsize=9,
        loc="upper right",
    )

    figure.suptitle(
        "Kaplan-Meier survival, five replicates pooled (n=1250 per curve)",
        color=INK,
        fontsize=12,
        x=0.006,
        ha="left",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output, dpi=DPI, facecolor=SURFACE)
    plt.close(figure)


def plot_replicates(per_replicate: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(
        1,
        len(POPULATIONS),
        figsize=(13, 4.0),
        sharey=True,
        sharex=True,
        facecolor=SURFACE,
    )

    for ax, population in zip(axes, POPULATIONS):
        apply_chrome(ax)
        colour = POPULATION_COLOURS[population]
        subset = per_replicate[per_replicate["population"] == population]

        for (_, _), curve in subset.groupby(["replicate", "sex"], sort=True):
            ax.step(
                curve.sort_values("age")["age"],
                curve.sort_values("age")["survival"],
                where="post",
                color=colour,
                linewidth=1.2,
                alpha=0.65,
                zorder=3,
            )

        ax.set_title(population, color=INK, fontsize=11, loc="left", pad=10)
        ax.set_xlabel("Age from eclosion (days)", color=INK_SECONDARY, fontsize=9)
        ax.set_ylim(0, 1.02)

    axes[0].set_xlim(0, per_replicate["age"].max() * 1.02)
    axes[0].set_ylabel("Surviving fraction", color=INK_SECONDARY, fontsize=9)
    figure.suptitle(
        "Survival by replicate: 5 replicates x 2 sexes per population, "
        "shared x-axis scale",
        color=INK,
        fontsize=12,
        x=0.006,
        ha="left",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(output, dpi=DPI, facecolor=SURFACE)
    plt.close(figure)


def validate_inputs(pooled: pd.DataFrame, per_replicate: pd.DataFrame) -> None:
    errors = []

    expected = set(POPULATIONS)
    for name, frame in (("pooled", pooled), ("per-replicate", per_replicate)):
        missing = expected - set(frame["population"].unique())
        if missing:
            errors.append(f"{name} curves missing population(s) {sorted(missing)}")

    unnamed = set(POPULATION_COLOURS) ^ expected
    if unnamed:
        errors.append(f"no colour assigned for population(s) {sorted(unnamed)}")

    if errors:
        raise ValidationError(
            "plotting failed validation:\n  - " + "\n  - ".join(errors)
        )


def main(results_dir: Path, out_dir: Path) -> None:
    pooled = pd.read_csv(results_dir / "survival_pooled.csv")
    per_replicate = pd.read_csv(results_dir / "survival_curves.csv")
    validate_inputs(pooled, per_replicate)

    out_dir.mkdir(parents=True, exist_ok=True)
    population_figure = out_dir / "survival_by_population.png"
    replicate_figure = out_dir / "survival_by_replicate.png"

    plot_populations(pooled, population_figure)
    plot_replicates(per_replicate, replicate_figure)

    print(f"wrote {population_figure}")
    print(f"wrote {replicate_figure}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "results_dir",
        type=Path,
        nargs="?",
        default=RESULTS_DIR,
        help="directory holding the survival curve tables from stage 4",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=RESULTS_DIR / "figures",
        help="directory for the PNG figures",
    )
    args = parser.parse_args()
    main(args.results_dir, args.out_dir)
