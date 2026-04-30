from __future__ import annotations

import pandas as pd

import _bootstrap  # noqa: F401
from config.config import (
    RESULTS_DIR,
    TCGA_SKCM_EXPRESSION_URL,
    TCGA_SKCM_PHENOTYPE_URL,
    TCGA_SKCM_SURVIVAL_URL,
)


def read_table(url: str) -> pd.DataFrame:
    return pd.read_csv(url, sep="\t")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rankings = pd.read_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t")
    top_tfs = rankings["tf"].head(20).tolist()

    expr = read_table(TCGA_SKCM_EXPRESSION_URL)
    phenotype = read_table(TCGA_SKCM_PHENOTYPE_URL)
    survival = read_table(TCGA_SKCM_SURVIVAL_URL)

    expr = expr.rename(columns={expr.columns[0]: "gene_id"})
    expr["gene_symbol"] = expr["gene_id"].astype(str).str.split("|").str[0]
    tf_expr = expr.loc[expr["gene_symbol"].isin(top_tfs)].copy()

    sample_columns = [col for col in tf_expr.columns if col not in {"gene_id", "gene_symbol"}]
    melted = tf_expr.melt(
        id_vars=["gene_id", "gene_symbol"],
        value_vars=sample_columns,
        var_name="sample",
        value_name="log2_fpkm_uq_plus_1",
    )
    summary = (
        melted.groupby("gene_symbol", as_index=False)["log2_fpkm_uq_plus_1"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )
    summary.columns = ["tf", "mean_bulk_expr", "median_bulk_expr", "std_bulk_expr", "n_samples"]

    phenotype.to_csv(RESULTS_DIR / "tcga_skcm_phenotype.tsv", sep="\t", index=False)
    survival.to_csv(RESULTS_DIR / "tcga_skcm_survival.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS_DIR / "tcga_skcm_tf_expression_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
