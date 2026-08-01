import os
from dataclasses import dataclass


@dataclass
class PipelinePaths:
    """Filesystem layout for the pipeline.

    Everything hangs off ``base``. ``scratch`` (large reference genomes) defaults
    to ``{base}/.scratch`` so it survives machine reboots -- do NOT point it at
    ``/tmp`` on hosts that wipe ``/tmp`` on restart. All fields can be overridden
    via constructor args or the environment variables below.
    """

    base: str = os.environ.get(
        "CATLAS_BASE", os.path.expanduser("~/hypothesis-generation_demo")
    )
    scratch: str = ""

    def __post_init__(self):
        if not self.scratch:
            self.scratch = os.environ.get("CATLAS_SCRATCH", f"{self.base}/.scratch")

    @property
    def reference_bim_prefix(self) -> str:

        return os.environ.get("CATLAS_BIM_PREFIX", f"{self.base}/data/EUR")

    @property
    def eur_bim_prefix(self) -> str:
        return self.reference_bim_prefix

    @property
    def fasta_hg38(self) -> str:
        return f"{self.scratch}/hg38.fa"

    @property
    def fasta_hg19(self) -> str:
        return f"{self.scratch}/hg19.fa"

    @property
    def liftover_bin(self) -> str:
        return os.environ.get("LIFTOVER_BIN", f"{self.base}/bin/liftOver")

    @property
    def chain_file(self) -> str:
        return f"{self.base}/data/reference/hg19ToHg38.over.chain.gz"

    @property
    def plink_bin(self) -> str:
        return os.environ.get("PLINK_BIN", os.path.expanduser("~/bin/plink"))

    @property
    def catlas_dir(self) -> str:
        return f"{self.base}/data/catlas_ccres"

    def dataset_dir(self, label: str) -> str:
        d = f"{self.base}/data/datasets/{label.lower()}"
        os.makedirs(d, exist_ok=True)
        return d

    def results_dir(self) -> str:
        d = f"{self.base}/data/results"
        os.makedirs(d, exist_ok=True)
        return d


@dataclass
class PipelineConfig:
    ld_prune_params: tuple = (
        "1000",
        "100",
        "0.1",
    )  # PLINK --indep-pairwise: window, step, r^2
    n_perm: int = 100_000
    perm_window_bp: int = 1_000_000
    promiscuity_threshold: float = 0.20  # exclude variants open in >20% of cell types
    fdr_threshold: float = 0.10
    seq_half_window_bp: int = 500  # +/-500bp, shared by NT-2.5B and Evo2
    build_check_n_test: int = 50  # variants sampled for build auto-detection
