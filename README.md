# sc-melanoma-celloracle

Melanoma scRNA-seq project scaffold for **GSE72056 / Tirosh et al. 2016** using **CellOracle** for GRN inference and in silico TF perturbation, with downstream validation against **TCGA-SKCM** bulk expression.

## Why CellOracle

- The available Tirosh matrix is **log2(TPM+1)**, not raw UMI counts.
- That makes it a poor fit for the official Geneformer tokenizer workflow.
- CellOracle is a better match because it operates on normalized scRNA-seq matrices and models TF perturbation through an inferred regulatory network.

## Current status

- Data staging and Scanpy QC run locally.
- The repo is now oriented around a **CellOracle-first** workflow.
- No valid CellOracle perturbation result has been committed yet.
- Earlier TF-IDF fallback outputs were removed from the scientific narrative because they were not Geneformer or CellOracle results.

## Dataset summary

- `GSE72056` metastatic melanoma single-cell RNA-seq from Tirosh et al. 2016.
- `4,645` total cells across `19` tumors.
- `1,257` malignant cells.
- `3,256` non-malignant cells.
- `132` unresolved cells.

## Workflow

1. Parse the Tirosh matrix into expression and metadata tables.
2. Build a filtered `AnnData` object with Scanpy.
3. Prepare a CellOracle-ready `AnnData` with embeddings and malignant-cell annotations.
4. Load a human promoter-based base GRN prior for CellOracle.
5. Infer regulatory links with CellOracle.
6. Simulate TF perturbations in malignant melanoma cells.
7. Validate shortlisted TFs against TCGA-SKCM bulk expression.

## Repo layout

- `config/config.py`: paths, URLs, thresholds, TF seed list, and CellOracle prior locations.
- `scripts/00_prepare_data.py`: download or stage GSE72056 and split metadata from expression.
- `scripts/01_qc_filter.py`: create and QC the main `AnnData`.
- `scripts/02_celloracle_preprocess.py`: prepare a CellOracle-ready `AnnData` with embeddings and malignant subsets.
- `scripts/03_celloracle_grn.py`: build the Oracle object, import the base GRN prior, and infer links.
- `scripts/04_celloracle_perturbation.py`: simulate TF perturbations and rank embedding shifts.
- `scripts/04_tcga_skcm_validation.py`: summarize TF expression in TCGA-SKCM bulk data.
- `scripts/05_build_report.py`: merge CellOracle perturbation outputs with TCGA summaries.

## Quick start

```bash
python scripts/00_prepare_data.py
python scripts/01_qc_filter.py
python scripts/02_celloracle_preprocess.py
python scripts/03_celloracle_grn.py --base-grn /path/to/celloracle_human_base_grn.parquet
python scripts/04_celloracle_perturbation.py --tf MITF --tf SOX10 --tf TFAP2A
python scripts/04_tcga_skcm_validation.py --source cbioportal \
  --expr refs/skcm_tcga_pan_can_atlas_2018_data_mrna_seq_v2_rsem.txt \
  --phenotype refs/skcm_tcga_pan_can_atlas_2018_data_clinical_sample.txt \
  --survival refs/skcm_tcga_pan_can_atlas_2018_data_clinical_patient.txt
python scripts/05_build_report.py
```

## Installation note

CellOracle installation on this local macOS arm64 environment currently fails in native dependencies such as `velocyto` and `gimmemotifs` because of compiler and OpenMP issues. The highest-confidence path is to run this workflow in a Linux or Docker-backed environment following the official CellOracle documentation.

## Scientific guardrails

- Do not describe this repo as having completed a Geneformer analysis on GSE72056.
- Do not interpret the earlier removed TF-IDF fallback outputs as biology.
- Treat TCGA validation as a downstream cross-check only after a real CellOracle perturbation result exists.

## Primary sources

- Tirosh et al. 2016, Science: [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056)
- CellOracle documentation: [official docs](https://morris-lab.github.io/CellOracle.documentation/)
- TCGA melanoma bulk cohort used here: cBioPortal SKCM PanCancer Atlas 2018 and related public clinical tables
