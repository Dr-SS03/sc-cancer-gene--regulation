from __future__ import annotations

import pandas as pd
import scanpy as sc
from anndata import AnnData

import _bootstrap  # noqa: F401
from config.config import MALIGNANT_LABEL, MAX_MT_PCT, MIN_CELLS, MIN_GENES, PROCESSED_DIR


def main() -> None:
    expr = pd.read_csv(PROCESSED_DIR / "expression_matrix.tsv", sep="\t", index_col=0).T
    meta = pd.read_csv(PROCESSED_DIR / "cell_metadata.tsv", sep="\t", index_col=0)

    adata = AnnData(X=expr.to_numpy(), obs=meta.copy(), var=pd.DataFrame(index=expr.columns))
    adata.obs_names = expr.index
    adata.var_names_make_unique()
    adata.obs["is_malignant"] = adata.obs["malignant"].astype(str).eq(MALIGNANT_LABEL)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    # Tirosh GSE72056 is already log2(TPM+1); QC here is a light filtering pass, not a raw-count UMI workflow.
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    adata = adata[adata.obs["pct_counts_mt"] < MAX_MT_PCT].copy()

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    adata.write(PROCESSED_DIR / "melanoma_qc.h5ad")
    print(adata)


if __name__ == "__main__":
    main()
