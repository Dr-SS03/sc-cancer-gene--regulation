from __future__ import annotations

import numpy as np
import scanpy as sc

import _bootstrap  # noqa: F401
from config.config import MALIGNANT_LABEL, PROCESSED_DIR


def main() -> None:
    adata = sc.read_h5ad(PROCESSED_DIR / "melanoma_qc.h5ad")
    adata.obs["malignant_label"] = adata.obs["malignant"].astype(str)
    adata.obs["celloracle_cluster"] = np.where(
        adata.obs["malignant_label"].eq(MALIGNANT_LABEL),
        adata.obs["tumor"].astype(str).radd("malignant_"),
        adata.obs["non_malignant_cell_type"].astype(str).radd("nonmal_"),
    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, resolution=0.6, key_added="celloracle_leiden")
    except ImportError:
        # Keep the preprocessing runnable even in lighter environments without igraph/leiden.
        adata.obs["celloracle_leiden"] = adata.obs["celloracle_cluster"].astype(str)

    adata.write(PROCESSED_DIR / "melanoma_celloracle_input.h5ad")
    adata[adata.obs["malignant_label"].eq(MALIGNANT_LABEL)].write(
        PROCESSED_DIR / "melanoma_malignant_celloracle_input.h5ad"
    )
    print(adata)


if __name__ == "__main__":
    main()
