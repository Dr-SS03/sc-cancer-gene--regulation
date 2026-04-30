from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import _bootstrap  # noqa: F401
from config.config import PROCESSED_DIR, RESULTS_DIR


def load_base_grn(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Base GRN prior not found: {path}")
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep=None, engine="python")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", default=str(PROCESSED_DIR / "melanoma_celloracle_input.h5ad"))
    parser.add_argument("--base-grn", required=True, help="Path to a CellOracle-compatible human base GRN prior.")
    args = parser.parse_args()

    try:
        import celloracle as co
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError(
            "CellOracle is not available in this environment. "
            "Use a Linux or Docker-backed environment that can install CellOracle cleanly."
        ) from exc

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.adata)
    base_grn = load_base_grn(Path(args.base_grn))

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="celloracle_cluster",
        embedding_name="X_umap",
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.perform_PCA()
    oracle.knn_imputation(
        n_pca_dims=min(30, adata.obsm["X_pca"].shape[1]),
        k=15,
        balanced=True,
        b_sight=120,
        b_maxl=60,
        n_jobs=1,
    )
    links = oracle.get_links(cluster_name_for_GRN_unit="celloracle_cluster", alpha=10, verbose_level=1)

    oracle.to_hdf5(file_path=str(RESULTS_DIR / "celloracle_oracle.hdf5"))
    links.to_hdf5(file_path=str(RESULTS_DIR / "celloracle_links.hdf5"))
    print("Saved CellOracle Oracle object and link table.")


if __name__ == "__main__":
    main()
