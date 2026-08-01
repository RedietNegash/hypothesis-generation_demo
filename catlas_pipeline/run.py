"""Command-line entry point for the CATlas cell-type enrichment pipeline.

Run a built-in dataset:
    python -m catlas_pipeline.run --preset t2d
    python -m catlas_pipeline.run --preset all          # t2d, bmi, af, ad

Run any custom GWAS by describing its columns (no code changes needed):
    python -m catlas_pipeline.run \
        --gwas /path/to/my_gwas.tsv.gz --label MYTRAIT \
        --col chr=CHR --col pos=BP --col ref=A1 --col alt=A2 \
        --col beta=BETA --col pval=P --col rsid=SNP \
        --pval 5e-8

If the file has no real dbSNP rsIDs (e.g. its id column is chr:pos), omit
--col rsid and add --backfill-rsids to recover them from the LD reference panel.
"""

import argparse

from . import ld
from .config import PipelineConfig, PipelinePaths
from .pipeline import ensure_reference_data, run_pipeline

PRESETS = ("t2d", "bmi", "af", "ad")


def prepare_preset(name: str, paths: PipelinePaths):
    """Download (if needed) and prepare one of the built-in datasets.

    These double as worked examples of column mapping + dataset-specific quirks.
    Returns (canonical_df, label).
    """
    from .catlas import download

    name = name.lower()
    if name == "bmi":
        from .prep.bmi import GWAS_URL, VAR_URL, prepare_bmi

        gwas = f"{paths.dataset_dir('BMI')}/gwas.tsv.bgz"
        var = f"{paths.dataset_dir('BMI')}/variants.tsv.bgz"
        download(GWAS_URL, gwas, "BMI GWAS (UK Biobank / Neale lab)")
        download(VAR_URL, var, "BMI variant annotation (rsID map)")
        return prepare_bmi(gwas, var), "BMI"

    if name == "af":
        from .prep.af import GWAS_URL, prepare_af

        gwas = f"{paths.dataset_dir('AF')}/gwas.tbl.gz"
        download(GWAS_URL, gwas, "AF GWAS (Nielsen/Thorolfsdottir/Willer 2018)")
        return prepare_af(gwas), "AF"

    if name == "t2d":
        from .prep.t2d import prepare_t2d

        gwas = f"{paths.dataset_dir('T2D')}/gwas.tsv.gz"
        df = prepare_t2d(gwas)
        # The Mahajan T2D file's "SNP" column is a chr:pos id, not a dbSNP rsID,
        # so it matches 0% of the reference panel -- recover real rsIDs by join.
        bim = ld.load_reference_bim(paths.reference_bim_prefix)
        df = ld.backfill_rsids_from_bim(df, bim, label="T2D")
        return df, "T2D"

    if name == "ad":
        from .prep.ad import GWAS_URL, prepare_ad

        gwas = f"{paths.dataset_dir('AD')}/gwas.txt.gz"
        download(GWAS_URL, gwas, "AD GWAS (Jansen et al. 2019)")
        return prepare_ad(gwas), "AD"

    raise ValueError(f"Unknown preset: {name!r} (choose from {PRESETS})")


def _parse_col(pairs):
    mapping = {}
    for p in pairs or []:
        if "=" not in p:
            raise ValueError(f"--col expects 'canonical=SOURCE', got {p!r}")
        key, val = p.split("=", 1)
        mapping[key.strip()] = val.strip()
    return mapping


def _report(label, results, shortlist, enrichment):
    print(f"\n\n########## [{label}] PIPELINE COMPLETE ##########")
    print(f"Shortlisted cell types ({len(shortlist)}): {shortlist}")
    print("\nTop 10 cell types by permutation p-value:")
    print(results.sort_values("p_perm").head(10).to_string(index=False))
    if enrichment is not None:
        print("\nSequence-model enrichment (in-peak vs out-of-peak):")
        print(enrichment.to_string(index=False))


def build_parser():
    ap = argparse.ArgumentParser(
        prog="python -m catlas_pipeline.run",
        description="CATlas cell-type enrichment pipeline for GWAS variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--preset", choices=PRESETS + ("all",),
        help="run a built-in dataset (downloads + prepares it)",
    )
    src.add_argument("--gwas", help="path to a custom GWAS summary-stats file")

    ap.add_argument("--label", help="dataset name for outputs (required with --gwas)")
    ap.add_argument(
        "--col", action="append", metavar="canonical=SOURCE",
        help="column mapping, repeatable. canonical in {chr,pos,ref,alt,pval,beta,"
             "rsid}. Required with --gwas: chr,pos,ref,alt,pval",
    )
    ap.add_argument("--sep", default="\t", help="column separator (default: TAB)")
    ap.add_argument("--pval", type=float, default=5e-8, help="p-value threshold")
    ap.add_argument("--uppercase-alleles", action="store_true",
                    help="upper-case the ref/alt allele strings")
    ap.add_argument(
        "--backfill-rsids", action="store_true",
        help="replace rsIDs via a chr:pos join against the reference panel "
             "(use when the file has no real dbSNP rsIDs)",
    )

    ap.add_argument("--base", help="project base dir "
                    "(default: $CATLAS_BASE or ~/hypothesis-generation_demo)")
    ap.add_argument("--scratch", help="scratch dir for reference genomes "
                    "(default: {base}/.scratch)")
    ap.add_argument("--no-scoring", action="store_true",
                    help="skip Evo2/NT-2.5B GPU scoring (screen + shortlist only)")
    ap.add_argument("--n-perm", type=int, help="permutations for the null "
                    "(default: 100000)")
    ap.add_argument("--fdr", type=float, help="FDR threshold for the shortlist "
                    "(default: 0.10)")
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)

    paths_kwargs = {}
    if args.base:
        paths_kwargs["base"] = args.base
    if args.scratch:
        paths_kwargs["scratch"] = args.scratch
    paths = PipelinePaths(**paths_kwargs)

    config = PipelineConfig()
    if args.n_perm:
        config.n_perm = args.n_perm
    if args.fdr is not None:
        config.fdr_threshold = args.fdr

    # Resolve the dataset job list up-front (validates args before heavy downloads).
    if args.preset == "all":
        jobs = [("preset", p) for p in PRESETS]
    elif args.preset:
        jobs = [("preset", args.preset)]
    else:
        if not args.label:
            ap.error("--label is required with --gwas")
        cmap = _parse_col(args.col)
        missing = [k for k in ("chr", "pos", "ref", "alt", "pval") if k not in cmap]
        if missing:
            ap.error(f"--col is missing required mappings: {missing}")
        jobs = [("custom", (args.gwas, args.label, cmap))]

    ensure_reference_data(paths)

    summary = []
    for kind, spec in jobs:
        if kind == "preset":
            df, label = prepare_preset(spec, paths)
        else:
            from .prep.generic import prepare_gwas

            gwas, label, cmap = spec
            df = prepare_gwas(
                gwas, cmap, pval_threshold=args.pval, sep=args.sep,
                uppercase_alleles=args.uppercase_alleles, label=label,
            )
            if args.backfill_rsids:
                bim = ld.load_reference_bim(paths.reference_bim_prefix)
                df = ld.backfill_rsids_from_bim(df, bim, label=label)

        results, shortlist, enrichment = run_pipeline(
            df, label=label, paths=paths, config=config,
            run_scoring=not args.no_scoring,
        )
        _report(label, results, shortlist, enrichment)
        summary.append((label, shortlist[:5] if shortlist else []))

    if len(summary) > 1:
        print("\n\n===== SUMMARY (top cell types per trait) =====")
        for label, top in summary:
            print(f"  {label}: {', '.join(top) if top else '(none)'}")


if __name__ == "__main__":
    main()
