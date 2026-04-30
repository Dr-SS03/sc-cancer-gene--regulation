# sc-melanoma-geneformer

Foundation-model-based scRNA-seq pipeline for **GSE72056 / Tirosh et al. 2016 melanoma** to identify transcription factors whose in silico perturbation shifts malignant melanoma cells away from a malignant state, then validate those TFs against **TCGA-SKCM** bulk expression.

## Project question

Can a Geneformer-based perturbation workflow prioritize transcription factors whose simulated knockout reduces the malignant transcriptional state of melanoma cells, and do those TFs also show supportive expression patterns in TCGA-SKCM bulk RNA-seq?

## Dataset choice

- `GSE72056` metastatic melanoma single-cell RNA-seq from Tirosh et al. 2016.
- Local source staged in this repo: `data/raw/GSE72056_melanoma_single_cell_revised_v2.txt.gz`
- Matrix format: genes on rows, cells on columns, with metadata embedded in the first three rows.
- Observed in the staged file:
- `4,645` total cells
- `19` tumors / patients
- `1,257` malignant cells (`malignant == 2`)
- `3,256` non-malignant cells (`malignant == 1`)
- `132` unresolved cells (`malignant == 0`)

## Analysis strategy

1. Parse the Tirosh matrix into a clean expression table and cell metadata.
2. Build an `AnnData` object and perform light QC / filtering with Scanpy.
3. Tokenize malignant and reference cells into Geneformer-compatible ranked gene lists.
4. Use Geneformer embeddings to define a malignant-state classifier.
5. Perform in silico TF perturbations across malignant cells and rank TFs by how strongly they reduce malignant-state probability.
6. Cross-check the highest-ranking TFs in TCGA-SKCM bulk RNA-seq from UCSC Xena / GDC-derived tables.

## Repo layout

- `config/config.py`: central paths, URLs, thresholds, and seed TF list.
- `scripts/00_prepare_data.py`: copy or download GSE72056, decompress it, split metadata and expression tables, and write a dataset summary.
- `scripts/01_qc_filter.py`: create a Scanpy `AnnData` object, annotate malignant status, and save a filtered `.h5ad`.
- `scripts/02_geneformer_inputs.py`: convert expression values into per-cell ranked gene lists and export malignant/reference metadata for Geneformer.
- `scripts/03_geneformer_perturbation.py`: compute Geneformer embeddings, fit a malignant-state classifier, and rank TF perturbations by predicted state shift.
- `scripts/04_tcga_skcm_validation.py`: download TCGA-SKCM expression / phenotype tables and summarize how top TFs behave in bulk tumors.
- `scripts/05_build_report.py`: merge perturbation and TCGA results into a publication-style candidate table.

## Quick start

```bash
conda env create -f environment.yml
conda activate sc-melanoma-geneformer

python scripts/00_prepare_data.py --local-gz data/raw/GSE72056_melanoma_single_cell_revised_v2.txt.gz
python scripts/01_qc_filter.py
python scripts/02_geneformer_inputs.py
python scripts/03_geneformer_perturbation.py
python scripts/04_tcga_skcm_validation.py
python scripts/05_build_report.py
```

## Important caveats

- The staged Tirosh matrix is a **preprocessed log2(TPM+1)** table, not raw UMI counts. Geneformer still works on ranked gene expression, but interpretation should stay framed as a ranking-based perturbation analysis rather than a raw-count reconstruction.
- The default Geneformer model on Hugging Face is newer than the original 2023 paper release. This repo pins Geneformer through the install command and keeps the workflow explicit in case you want to swap to a different released checkpoint later.
- `scripts/03_geneformer_perturbation.py` assumes Geneformer and Torch are installed and that model weights can be fetched. CPU execution is possible for this dataset, but it will be slower than GPU.
- TCGA-SKCM validation in this repo is intentionally lightweight and transparent: it checks whether top TFs are expressed in melanoma bulk tumors and whether expression differs by clinically available groups when those columns are present.

## Primary sources

- Tirosh et al. 2016, Science: [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056)
- Geneformer model card and docs: [Hugging Face](https://huggingface.co/ctheodoris/Geneformer), [documentation](https://geneformer.readthedocs.io/en/latest/)
- TCGA / UCSC Xena: [UCSC Xena public data](https://xena.ucsc.edu/public/), [GDC data types](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/types)

## Status

This repo has been scaffolded around the exact melanoma-to-Geneformer-to-TCGA question. The data parsing and project structure are ready; the Geneformer and TCGA scripts are prepared for execution once the full environment is installed.
