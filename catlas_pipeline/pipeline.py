import os

import pandas as pd

from . import enrichment, genome, ld, scoring
from .catlas import (
    CELL_TYPES,
    download_catlas_celltypes,
    filter_promiscuous_variants,
    get_shortlist,
    screen_celltypes,
)
from .config import PipelineConfig, PipelinePaths


def ensure_reference_data(paths: PipelinePaths, cell_types=None):

    os.makedirs(paths.scratch, exist_ok=True)
    os.makedirs(f"{paths.base}/data/reference", exist_ok=True)
    genome.ensure_reference_genomes(paths)
    ready_cell_types = download_catlas_celltypes(
        paths.catlas_dir, cell_types or CELL_TYPES
    )
    return ready_cell_types


def run_pipeline(
    df: pd.DataFrame,
    label: str,
    paths: PipelinePaths = None,
    config: PipelineConfig = None,
    run_scoring: bool = True,
    cell_types: list = None,
):

    paths = paths or PipelinePaths()
    config = config or PipelineConfig()
    cell_types = cell_types or CELL_TYPES
    dataset_dir = paths.dataset_dir(label)

    print(f"\n=== [{label}] Step 1: {len(df):,} variants loaded ===")

    print(f"\n=== [{label}] Step 2: rsID matching ===")
    bim = ld.load_eur_bim(paths.eur_bim_prefix)
    matched_df = ld.match_rsids(df, set(bim["rsid"]), label=label)

    print(f"\n=== [{label}] Step 3: LD pruning ===")
    pruned_df = ld.ld_prune(
        matched_df,
        paths.eur_bim_prefix,
        paths.plink_bin,
        dataset_dir,
        config.ld_prune_params,
        label=label,
    )

    print(f"\n=== [{label}] Step 3.5: genome build ===")
    genome_hg19 = genome.load_genome(paths.fasta_hg19)
    genome_hg38 = genome.load_genome(paths.fasta_hg38)
    build = genome.detect_genome_build(
        pruned_df, genome_hg19, genome_hg38, n_test=config.build_check_n_test
    )
    print(f"[{label}] Detected build: {build}")

    if build == "hg19":
        df_hg38 = genome.liftover_hg19_to_hg38(
            pruned_df,
            paths.liftover_bin,
            paths.chain_file,
            dataset_dir,
            tag=label.lower(),
        )
    else:
        df_hg38 = pruned_df.copy()

    mismatches = genome.verify_against_hg38(
        df_hg38, genome_hg38, label=f"{label} post-liftover"
    )
    assert mismatches == 0, (
        f"[{label}] Liftover verification failed ({mismatches} mismatches) -- "
        f"do not proceed to sequence extraction"
    )

    df_hg38 = genome.reorient_alleles(df_hg38, genome_hg38)

    print(f"\n=== [{label}] Promiscuity filter ===")
    df_filtered = filter_promiscuous_variants(
        df_hg38,
        paths.catlas_dir,
        config.promiscuity_threshold,
        label=label,
        cell_types=cell_types,
    )

    print(f"\n=== [{label}] Step 4: sequence extraction ===")
    seq_path = f"{dataset_dir}/sequences.tsv"
    seq_df = genome.extract_sequences(
        df_filtered, genome_hg38, config.seq_half_window_bp, seq_path
    )

    nt_scores, evo2_scores, nt_embeddings = None, None, None
    if run_scoring:
        print(f"\n=== [{label}] Step 5a: Evo2 scoring ===")
        evo2_scores = scoring.score_with_evo2(seq_df, f"{dataset_dir}/scores_evo2.tsv")

        print(f"\n=== [{label}] Step 5b: NT-2.5B scoring ===")
        nt_scores = scoring.score_with_nt(
            seq_df, f"{dataset_dir}/scores_nt.tsv", f"{dataset_dir}/nt_embeddings.npy"
        )
        import numpy as np

        nt_embeddings = np.load(f"{dataset_dir}/nt_embeddings.npy")

    print(f"\n=== [{label}] Step 6: CATlas cell-type screen ===")
    results_222 = screen_celltypes(
        df_filtered,
        paths.catlas_dir,
        label=label,
        n_perm=config.n_perm,
        window=config.perm_window_bp,
        cell_types=cell_types,
    )
    results_222.to_csv(
        f"{paths.results_dir()}/celltype_screen_222_{label.lower()}.tsv",
        sep="\t",
        index=False,
    )

    shortlist = get_shortlist(results_222, fdr_threshold=config.fdr_threshold)

    enrichment_df = None
    if run_scoring and shortlist:
        print(f"\n=== [{label}] Step 7: in-peak vs out-of-peak ===")
        enrichment_df = enrichment.compare_in_out_peak(
            nt_scores,
            evo2_scores,
            nt_embeddings,
            shortlist,
            results_222,
            paths.catlas_dir,
            f"{paths.results_dir()}/{label.lower()}_celltype_enrichment.tsv",
            label=label,
        )

    return results_222, shortlist, enrichment_df
