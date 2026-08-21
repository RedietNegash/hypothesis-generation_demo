import math
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.power import TTestIndPower


BASE_DIR = Path("/mnt/hdd_1/rediet/fly-ldsc")
GWAS_FILE = BASE_DIR / "data" / "gwas" / "lifespan_female_gwas.tsv"
FREQ_FILE = BASE_DIR / "data" / "gwas" / "tmp" / "female_lifespan_freq.frq"
POWER_FILE = BASE_DIR / "data" / "gwas" / "tmp" / "female_lifespan_power.tsv"
POWER_SUMMARY_FILE = BASE_DIR / "data" / "gwas" / "tmp" / "female_lifespan_power_summary.tsv"

N = 197
SD = 9.895947797640865
ALPHA = 2.28e-8

EFFECTS = [1.0, 2.0, 5.0, 7.5, 10.0, 15.0, 20.0]
MAFS = np.arange(0.01, 0.51, 0.01)


def homozygote_power(n, maf, effect_days, sd, alpha):
    n_minor_hom = n * maf**2
    n_common_hom = n * (1 - maf)**2

    if n_minor_hom < 1 or n_common_hom < 1:
        return np.nan

    pooled_sd = sd
    d = abs(effect_days) / pooled_sd

    analysis = TTestIndPower()

    return analysis.solve_power(
        effect_size=d,
        nobs1=n_common_hom,
        alpha=alpha,
        power=None,
        ratio=n_minor_hom / n_common_hom,
        alternative="two-sided",
    )


def run_power_analysis():
    rows = []

    for maf in MAFS:
        for effect in EFFECTS:
            power = homozygote_power(
                N,
                maf,
                effect,
                SD,
                ALPHA,
            )

            rows.append({
                "N": N,
                "MAF": round(maf, 2),
                "effect_days": effect,
                "sd_lifespan": SD,
                "alpha": ALPHA,
                "n_common_homozygotes": N * (1 - maf)**2,
                "n_minor_homozygotes": N * maf**2,
                "power": power,
            })

    df = pd.DataFrame(rows)

    POWER_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(POWER_FILE, sep="\t", index=False)

    print(f"Power results: {POWER_FILE}")

    summary = []

    for effect in EFFECTS:
        subset = df[df["effect_days"] == effect].dropna(subset=["power"])

        eligible = subset[subset["power"] >= 0.80]

        if len(eligible) > 0:
            min_maf = eligible["MAF"].min()
        else:
            min_maf = np.nan

        max_power = subset["power"].max()

        summary.append({
            "effect_days": effect,
            "minimum_MAF_for_80pct_power": min_maf,
            "maximum_power": max_power,
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(
        POWER_SUMMARY_FILE,
        sep="\t",
        index=False,
    )

    print(f"Power summary: {POWER_SUMMARY_FILE}")


if __name__ == "__main__":
    run_power_analysis()