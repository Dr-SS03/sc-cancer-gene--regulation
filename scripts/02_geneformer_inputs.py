from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import _bootstrap  # noqa: F401
from config.config import MALIGNANT_LABEL, PROCESSED_DIR, REFERENCE_LABEL, RESULTS_DIR


def rank_genes_for_cell(values: np.ndarray, genes: np.ndarray, top_n: int = 2048) -> list[str]:
    nonzero = values > 0
    if not np.any(nonzero):
        return []
    ranked_idx = np.argsort(values[nonzero])[::-1][:top_n]
    return genes[nonzero][ranked_idx].tolist()


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROCESSED_DIR / "melanoma_qc.h5ad")
    keep = adata.obs["malignant"].astype(str).isin([MALIGNANT_LABEL, REFERENCE_LABEL])
    adata = adata[keep].copy()

    rows = []
    genes = adata.var_names.to_numpy()
    for idx, cell_id in enumerate(adata.obs_names):
        ranked_genes = rank_genes_for_cell(np.asarray(adata.X[idx]).ravel(), genes)
        rows.append(
            {
                "cell_id": cell_id,
                "tumor": adata.obs.iloc[idx]["tumor"],
                "malignant": str(adata.obs.iloc[idx]["malignant"]),
                "ranked_genes": " ".join(ranked_genes),
            }
        )

    out = pd.DataFrame(rows)
    out.to_parquet(RESULTS_DIR / "geneformer_ranked_cells.parquet", index=False)
    out[["cell_id", "tumor", "malignant"]].to_csv(
        RESULTS_DIR / "geneformer_cell_metadata.tsv", sep="\t", index=False
    )
    print(f"Wrote {len(out)} ranked transcriptomes.")


if __name__ == "__main__":
    main()
