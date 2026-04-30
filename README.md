# sc-melanoma-geneformer

Foundation-model-based scRNA-seq pipeline for **GSE72056 / Tirosh et al. 2016 melanoma** to identify transcription factors whose in silico perturbation shifts malignant melanoma cells away from a malignant state, then validate those TFs against **TCGA-SKCM** bulk expression.

## Project question

Can a Geneformer-based perturbation workflow prioritize transcription factors whose simulated knockout reduces the malignant transcriptional state of melanoma cells, and do those TFs also show supportive expression patterns in TCGA-SKCM bulk RNA-seq?

## Dataset choice

- `GSE72056` metastatic melanoma single-cell RNA-seq from Tirosh et al. 2016.
- Local source used during development: `GSE72056_melanoma_single_cell_revised_v2.txt.gz`
- The raw GEO download is **not tracked in Git**. Use `scripts/00_prepare_data.py` to download it directly or point the script to a local copy.
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
3. Convert malignant and reference cells into ranked gene programs.
4. Run one of two Geneformer-compatible routes:
5. `official`: requires raw-count `.h5ad` or `.loom` with `ensembl_id` and `n_counts`, matching the current Geneformer tokenizer documentation.
6. `ranked_fallback`: uses the Tirosh log2(TPM+1) matrix as a rank-based approximation and scores TF deletions with a transparent surrogate embedding model.
7. Rank TF perturbations by how strongly they reduce malignant-state probability.
8. Cross-check the highest-ranking TFs in TCGA-SKCM bulk RNA-seq from UCSC Xena / GDC-derived tables.

## Repo layout

- `config/config.py`: central paths, URLs, thresholds, and seed TF list.
- `scripts/00_prepare_data.py`: copy or download GSE72056, decompress it, split metadata and expression tables, and write a dataset summary.
- `scripts/01_qc_filter.py`: create a Scanpy `AnnData` object, annotate malignant status, and save a filtered `.h5ad`.
- `scripts/02_geneformer_inputs.py`: convert expression values into per-cell ranked gene lists and export malignant/reference metadata.
- `scripts/03_geneformer_perturbation.py`: run an official Geneformer path when valid tokenizer inputs are available, or a rank-based fallback on the staged Tirosh matrix.
- `scripts/04_tcga_skcm_validation.py`: download TCGA-SKCM expression / phenotype tables, summarize top TF bulk behavior, and write plotting tables.
- `scripts/05_build_report.py`: merge perturbation and TCGA results into a publication-style candidate table.

## Quick start

```bash
conda env create -f environment.yml
conda activate sc-melanoma-geneformer

python scripts/00_prepare_data.py
python scripts/01_qc_filter.py
python scripts/02_geneformer_inputs.py
python scripts/03_geneformer_perturbation.py --mode ranked_fallback
python scripts/04_tcga_skcm_validation.py
python scripts/05_build_report.py
```

## Important caveats

- The staged Tirosh matrix is a **preprocessed log2(TPM+1)** table, not raw UMI counts. According to the current Geneformer tokenizer documentation, the official tokenizer expects raw-count `.h5ad` or `.loom` with `ensembl_id` and `n_counts`.
- For that reason, this repo exposes two perturbation modes. `official` is the standards-compliant Geneformer route if you later supply raw counts; `ranked_fallback` is the practical route for the exact Tirosh matrix you downloaded here.
- The default Geneformer model on Hugging Face is newer than the original 2023 paper release. This repo pins Geneformer through the install command and keeps the workflow explicit in case you want to swap to a different released checkpoint later.
- TCGA-SKCM validation in this repo now produces sample-level and group-level tables designed for downstream plots, but it still depends on the columns exposed by the UCSC Xena phenotype release you download at run time.

## Primary sources

- Tirosh et al. 2016, Science: [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056)
- Geneformer model card and docs: [Hugging Face](https://huggingface.co/ctheodoris/Geneformer), [documentation](https://geneformer.readthedocs.io/en/latest/)
- TCGA / UCSC Xena: [UCSC Xena public data](https://xena.ucsc.edu/public/), [GDC data types](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/types)

## Status

This repo is now structured as a cleaner public project: raw GEO input files stay outside Git history, the Geneformer step explicitly distinguishes official versus fallback workflows, and the TCGA validation step is set up to produce analysis-ready summaries.
