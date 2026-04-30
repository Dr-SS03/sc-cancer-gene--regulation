from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
REFS_DIR = PROJECT_ROOT / "refs"

GSE72056_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/"
    "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
)

TCGA_SKCM_EXPRESSION_URL = (
    "https://gdc.xenahubs.net/download/TCGA-SKCM.htseq_fpkm-uq.tsv.gz"
)
TCGA_SKCM_PHENOTYPE_URL = (
    "https://gdc.xenahubs.net/download/TCGA-SKCM.GDC_phenotype.tsv.gz"
)
TCGA_SKCM_SURVIVAL_URL = (
    "https://gdc.xenahubs.net/download/TCGA-SKCM.survival.tsv"
)
TCGA_SKCM_EXPRESSION_LOCAL = REFS_DIR / "TCGA-SKCM.htseq_fpkm-uq.tsv.gz"
TCGA_SKCM_PHENOTYPE_LOCAL = REFS_DIR / "TCGA-SKCM.GDC_phenotype.tsv.gz"
TCGA_SKCM_SURVIVAL_LOCAL = REFS_DIR / "TCGA-SKCM.survival.tsv"
TCGA_SKCM_SAMPLE_TYPE_CANDIDATES = [
    "_sample_type",
    "sample_type",
    "sample_type.samples",
    "definition",
]
TCGA_SKCM_PRIMARY_METASTATIC_CANDIDATES = [
    "_sample_type",
    "sample_type",
    "sample_type.samples",
]
CELLORACLE_BASE_GRN_LOCAL = REFS_DIR / "celloracle_human_base_grn.parquet"
CELLORACLE_DOCKER_NOTE = (
    "CellOracle native installation is not currently working in this local macOS arm64 environment. "
    "Prefer Linux or Docker for the actual GRN and perturbation run."
)

MIN_GENES = 200
MIN_CELLS = 3
MAX_MT_PCT = 20.0
MALIGNANT_LABEL = "2"
REFERENCE_LABEL = "1"
UNRESOLVED_LABEL = "0"

SEED_TFS = [
    "MITF",
    "SOX10",
    "TFAP2A",
    "TFAP2C",
    "STAT1",
    "STAT3",
    "IRF1",
    "IRF4",
    "JUN",
    "JUNB",
    "FOS",
    "RELA",
    "NFKB1",
    "E2F1",
    "MYC",
    "FOXM1",
    "TEAD1",
    "TEAD4",
    "YAP1",
    "AHR",
    "SOX9",
    "ATF4",
    "KLF4",
    "SMAD3",
]
