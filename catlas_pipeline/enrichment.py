import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .catlas import build_variant_pyranges, get_celltype_peaks


def compare_in_out_peak(
    nt_scores: pd.DataFrame,
    evo2_scores: pd.DataFrame,
    nt_embeddings: np.ndarray,
    shortlist: list,
    results_222: pd.DataFrame,
    catlas_dir: str,
    out_path: str,
    label: str = "",
) -> pd.DataFrame:
    merged = nt_scores.copy()
    if evo2_scores is not None:
        merged = merged.merge(
            evo2_scores.drop(columns=["pval"], errors="ignore"),
            on=["variant", "chr", "pos"],
            how="outer",
        )
    merged = merged.reset_index(drop=True)
    print(f"[{label}] Merged variant scores: {len(merged):,} rows")

    merged_pr = build_variant_pyranges(merged, keep_idx=True)

    sig_lookup = results_222.set_index("cell_type")[["p_perm", "FDR"]]

    results = []
    for cell_type in shortlist:
        peak_pr = get_celltype_peaks(cell_type, catlas_dir)
        if peak_pr is None:
            continue
        ov = merged_pr.overlap(peak_pr)
        in_idx = set(ov.variant_idx) if len(ov) else set()
        in_peak = merged.index.isin(in_idx)
        n_in, n_total = int(in_peak.sum()), len(merged)
        if n_in == 0:
            print(f"{cell_type}: 0 variants overlap peaks -- skipping")
            continue

        if (
            in_peak.sum()
            and (~in_peak).sum()
            and nt_embeddings is not None
            and nt_embeddings.shape[0] == len(merged)
        ):
            from scipy.spatial.distance import cosine

            nt_enrich = float(
                cosine(
                    nt_embeddings[in_peak].mean(axis=0),
                    nt_embeddings[~in_peak].mean(axis=0),
                )
            )
        else:
            nt_enrich = np.nan

        if "evo2_delta" in merged.columns:
            in_evo = merged.loc[in_peak, "evo2_delta"].abs().dropna()
            out_evo = merged.loc[~in_peak, "evo2_delta"].abs().dropna()
            if len(in_evo) and len(out_evo):
                _, p_evo2 = mannwhitneyu(in_evo, out_evo, alternative="greater")
                evo2_delta = in_evo.mean() - out_evo.mean()
            else:
                p_evo2, evo2_delta = np.nan, np.nan
        else:
            p_evo2, evo2_delta = np.nan, np.nan

        results.append(
            {
                "cell_type": cell_type,
                "n_variants_in_peaks": n_in,
                "n_total_variants": n_total,
                "p_perm": sig_lookup.loc[cell_type, "p_perm"],
                "FDR": sig_lookup.loc[cell_type, "FDR"],
                "nt_embedding_cosine_dist": round(nt_enrich, 6)
                if not np.isnan(nt_enrich)
                else np.nan,
                "evo2_delta": round(evo2_delta, 6)
                if not np.isnan(evo2_delta)
                else np.nan,
                "p_evo2": p_evo2,
            }
        )
        print(f"{cell_type}: {n_in:,}/{n_total:,} variants in accessible peaks")

    results_df = (
        pd.DataFrame(results).sort_values("p_perm") if results else pd.DataFrame()
    )
    if len(results_df):
        results_df["fdr_significant"] = results_df["FDR"] < 0.10
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved -> {out_path}")
    if len(results_df) and not results_df["fdr_significant"].any():
        print(
            "NOTE: no cell type in this table passes FDR<0.10 -- this is an "
            "exploratory ranking, not a confirmed hit list."
        )
    print(results_df.to_string(index=False))
    return results_df
