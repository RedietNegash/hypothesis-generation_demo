# Fly lifespan pipeline

Calculates lifespan data from the raw daily death counts in
`counts-8nqxycmkc7ysxro9d1ucfisuio.xlsx`, for the three experimentally evolved
_Drosophila melanogaster_ populations described in `superFlyProgressReport.html`:
**B** (2-week generation), **O** (10-week) and **SO**, super-O (16-week), five
replicate populations each, 250 males and 250 females per replicate.

## Running it

```bash
pip install -r requirements.txt
python run_all.py                    # all six stages, default workbook
python run_all.py path/to/other.xlsx # a different workbook
```

Every stage also runs standalone and reads the previous stage's default output,
so a single stage can be re-run after an edit without redoing the workbook
parse. Outputs land in `data/` and `results/`, both resolved against this
directory rather than the working directory.

## Stages

| #   | Script                    | Produces                                                                                          |
| --- | ------------------------- | ------------------------------------------------------------------------------------------------- |
| 1   | `load_counts.py`          | `data/deaths_long.csv` — one row per population/replicate/sex/day                                 |
| 2   | `expand_ages.py`          | `data/flies.csv` — one row per fly with its age at death, 7,500 rows                              |
| 3   | `summarise_replicates.py` | `lifespan_summary.csv`, `replicate_summary.csv`, `replicate_spread.csv`                           |
| 4   | `survival_curves.py`      | `survival_pooled.csv`, `survival_curves.csv`, `mortality.csv`, `mortality_replicate.csv`          |
| 5   | `compare_populations.py`  | `population_comparison.csv`, `sex_comparison.csv`, `replicate_measures.csv`, `pooled_logrank.csv` |
| 6   | `plot_survival.py`        | `figures/survival_by_population.png`, `figures/survival_by_replicate.png`                         |

`results/lifespan_summary.csv` is the headline table. Everything else is either
finer-grained or a curve.

Each stage validates its own output and raises `ValidationError` rather than
writing a table that a later stage would silently trust.

## The lifespan data

Age is counted in days from eclosion. All figures pool the five replicates.

| Population | Sex     | n    | Mean   | SD    | Median | KM median | Oldest | Top 10% mean |
| ---------- | ------- | ---- | ------ | ----- | ------ | --------- | ------ | ------------ |
| B          | Females | 1250 | 31.18  | 14.39 | 32     | 32        | 61     | 53.1         |
| B          | Males   | 1250 | 30.83  | 14.14 | 31     | 31        | 60     | 52.7         |
| O          | Females | 1250 | 70.24  | 18.69 | 72     | 72        | 102    | 96.7         |
| O          | Males   | 1250 | 69.18  | 19.38 | 72     | 72        | 102    | 96.4         |
| SO         | Females | 1250 | 135.23 | 51.27 | 141    | 141       | 225    | 211.2        |
| SO         | Males   | 1250 | 134.83 | 49.12 | 139    | 139       | 225    | 210.2        |

`top10_mean`, the mean of the longest-lived 10%, stands in for maximum lifespan;
the single oldest fly is one draw and moves with luck.

## Comparing B, O and SO

The five replicate populations are the experimental units, not the 250 flies
inside each. Flies in one replicate share a vial, a food batch and a handling
schedule, so their ages at death are correlated, and a test over pooled flies
reports a precision the design never bought. Stage 5 reduces each replicate to
one number and tests across the five with Welch's t-test, so n = 5 per group.

| Pair      | Ratio | Difference (days) | p (females) | p (males) |
| --------- | ----- | ----------------- | ----------- | --------- |
| O over B  | 2.25× | +39.1             | 1.2×10⁻¹⁰   | 4.6×10⁻¹⁰ |
| SO over O | 1.93× | +65.0             | 3.4×10⁻¹³   | 3.7×10⁻¹⁰ |
| SO over B | 4.34× | +104.1            | 3.7×10⁻¹⁴   | 1.2×10⁻¹¹ |

The same ordering holds for median lifespan and for the top-10% mean, in both
sexes. Replicates agree closely — the coefficient of variation of the five
replicate means is under 6% in every group and lowest in SO.

**Sex makes no detectable difference.** Females outlive males by 0.35 to 1.07
days, paired within replicate, p = 0.34 to 0.77 in all three populations.

`pooled_logrank.csv` runs the same pairs over the 1,250 pooled flies and is
recorded for contrast only: it returns chi-squared near 2,200 with p underflowing
to zero, against p ≈ 10⁻¹⁰ from the replicate-level test. Same conclusion, since
the effect is enormous, but the pooled figures are not quotable.

## Figures

`results/figures/survival_by_population.png` — Kaplan-Meier curves for the three
populations with 95% confidence bands, one panel per sex, five replicates pooled.

`results/figures/survival_by_replicate.png` — all 10 curves per population
(5 replicates x 2 sexes) on a shared x-axis, so both the tight agreement within a
population and the 4.5x spread between populations are visible at once.

## Age origin, and reading the progress report

`Day 1` is 2019-05-13 in all three sheets and the columns are already split by
sex, so the flies were sexed adults when counting began: **Day 1 is eclosion,
not egg lay**. `age_at_death` is therefore adult age, the standard for a
Drosophila lifespan assay.

The progress report quotes average lifespans of **40 / 80 / 140 days** for
B / O / SO. This data gives **31.0 / 69.7 / 135.0**. The gap is 9-10 days in
every population, which is the egg-to-eclosion development time at 25 °C: the
report counts from egg, and the two agree once that is added. `AGE_OFFSET_DAYS`
in `load_counts.py` shifts every reported age if you want the from-egg
convention; it defaults to 0.

## The workbook

Three `<population> death` sheets, each `Date, Day`, then five `Males`/`Females`
pairs. Checked against the data, not assumed:

- Every one of the 30 replicate×sex cohorts sums to **exactly 250** deaths.
- **Nothing is censored.** All 250 flies died within the recorded span, so every
  fly carries `event_observed = 1` and survival reaches zero. SO's last row
  (day 225) still has deaths on it, which looks like truncation, but the totals
  confirm the cohorts were exhausted.
- Blanks appear **only after** a cohort's last death, where they mean zero. A
  blank between recorded deaths would be a missed census whose flies died later,
  so stage 1 rejects the workbook rather than reading such a blank as zero.
- `Day` runs consecutively from 1 and the calendar dates agree with it.

The `<population> eggs` sheets hold weekly fecundity counts. This pipeline does
not read them.

## Not covered

- **No aging-rate model.** Nothing fits a Gompertz curve to the hazard, so there
  is no aging rate or mortality-rate doubling time.
- **No fecundity analysis.** The `eggs` sheets are untouched.
