# sc-melanoma-geneformer

Reproducible melanoma scRNA-seq project scaffold for **GSE72056 / Tirosh et al. 2016** with explicit preparation steps for a future Geneformer perturbation analysis and a bulk-expression validation layer against **TCGA-SKCM**.

## Project question

Can a Geneformer-based perturbation workflow prioritize transcription factors whose simulated knockout reduces the malignant transcriptional state of melanoma cells, and do those TFs also show supportive expression patterns in TCGA-SKCM bulk RNA-seq?

At the moment, this repo does **not** answer that question with a valid Geneformer run on GSE72056, because the available Tirosh matrix is log2(TPM+1) rather than tokenizer-compliant raw counts.

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
4. Prepare an `official` Geneformer route that requires raw-count `.h5ad` or `.loom` with `ensembl_id` and `n_counts`, matching the current Geneformer tokenizer documentation.
5. Provide a clearly labeled `ranked_fallback` diagnostic baseline for software plumbing only.
6. Cross-check candidate TFs in TCGA-SKCM bulk RNA-seq once a biologically valid perturbation signal exists.

## Repo layout

- `config/config.py`: central paths, URLs, thresholds, and seed TF list.
- `scripts/00_prepare_data.py`: copy or download GSE72056, decompress it, split metadata and expression tables, and write a dataset summary.
- `scripts/01_qc_filter.py`: create a Scanpy `AnnData` object, annotate malignant status, and save a filtered `.h5ad`.
- `scripts/02_geneformer_inputs.py`: convert expression values into per-cell ranked gene lists and export malignant/reference metadata.
- `scripts/03_geneformer_perturbation.py`: contains the unimplemented official Geneformer entry point plus a diagnostic fallback baseline that should not be interpreted biologically.
- `scripts/04_tcga_skcm_validation.py`: download TCGA-SKCM expression / phenotype tables, summarize TF bulk behavior, and write plotting tables.
- `scripts/05_build_report.py`: merge perturbation and TCGA summaries into a diagnostic table and fail fast on unsupported modes.

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
- For that reason, the executed `ranked_fallback` mode is **not Geneformer** and should not be described as a foundation-model perturbation analysis.
- The `ranked_fallback` output is a TF-IDF plus logistic-regression baseline built from ranked gene tokens. Removing one TF token from a 2048-token document produces extremely small score shifts, so those perturbation values are best treated as software diagnostics rather than biology.
- Any ranking that multiplies fallback perturbation scores by melanoma bulk expression is misleading and should not be used for candidate nomination.
- The default Geneformer model on Hugging Face is newer than the original 2023 paper release. This repo pins Geneformer through the install command and keeps the workflow explicit in case you want to swap to a different released checkpoint later.
- TCGA-SKCM validation in this repo now produces sample-level and group-level tables designed for downstream plots, but it still depends on the columns exposed by the UCSC Xena phenotype release you download at run time.

## Primary sources

- Tirosh et al. 2016, Science: [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056)
- Geneformer model card and docs: [Hugging Face](https://huggingface.co/ctheodoris/Geneformer), [documentation](https://geneformer.readthedocs.io/en/latest/)
- TCGA / UCSC Xena: [UCSC Xena public data](https://xena.ucsc.edu/public/), [GDC data types](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/types)

## Status

This repo is structured as a cleaner public scaffold: raw GEO input files stay outside Git history, the official Geneformer path is explicitly marked as pending valid raw-count inputs, and the TCGA validation layer is ready once a real perturbation result exists.
