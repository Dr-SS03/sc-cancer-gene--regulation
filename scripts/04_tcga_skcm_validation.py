from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import _bootstrap  # noqa: F401
from config.config import (
    FIGURES_DIR,
    RESULTS_DIR,
    REFS_DIR,
    TCGA_SKCM_PRIMARY_METASTATIC_CANDIDATES,
    TCGA_SKCM_EXPRESSION_URL,
    TCGA_SKCM_EXPRESSION_LOCAL,
    TCGA_SKCM_PHENOTYPE_LOCAL,
    TCGA_SKCM_PHENOTYPE_URL,
    TCGA_SKCM_SAMPLE_TYPE_CANDIDATES,
    TCGA_SKCM_SURVIVAL_LOCAL,
    TCGA_SKCM_SURVIVAL_URL,
)


def read_table(path_or_url: str | Path) -> pd.DataFrame:
    return pd.read_csv(path_or_url, sep="\t")


def read_cbioportal_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def find_first_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in frame.columns:
            return col
    return None


def classify_cbioportal_sample_type(value: object) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).lower()
    if "metast" in text:
        return "Metastatic"
    if "primary" in text:
        return "Primary"
    return None


def classify_primary_vs_metastatic(value: object) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).lower()
    if "metastatic" in text:
        return "Metastatic"
    if "primary" in text:
        return "Primary"
    return None


def plot_top_tfs(summary: pd.DataFrame, out_path: Path) -> None:
    plot_df = summary.sort_values("mean_bulk_expr", ascending=False).head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(plot_df["tf"], plot_df["mean_bulk_expr"], color="#2b6cb0")
    ax.set_xlabel("Mean TCGA-SKCM log2(FPKM-UQ + 1)")
    ax.set_ylabel("Transcription factor")
    ax.set_title("Top perturbation TFs in TCGA-SKCM bulk")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_primary_vs_metastatic(group_summary: pd.DataFrame, out_path: Path) -> None:
    if group_summary.empty:
        return
    pivot = group_summary.pivot(index="tf", columns="disease_group", values="mean_log2_fpkm_uq_plus_1")
    pivot = pivot.dropna()
    if pivot.empty:
        return
    pivot = pivot.sort_values("Metastatic", ascending=True).tail(15)
    fig, ax = plt.subplots(figsize=(8, 6))
    width = 0.38
    y = range(len(pivot))
    ax.barh([v - width / 2 for v in y], pivot["Primary"], height=width, label="Primary", color="#68a357")
    ax.barh([v + width / 2 for v in y], pivot["Metastatic"], height=width, label="Metastatic", color="#c05621")
    ax.set_yticks(list(y))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("Mean TCGA-SKCM log2(FPKM-UQ + 1)")
    ax.set_ylabel("Transcription factor")
    ax.set_title("Primary vs metastatic TCGA-SKCM expression")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of TFs from the perturbation ranking to carry into TCGA validation.",
    )
    parser.add_argument("--expr", default=str(TCGA_SKCM_EXPRESSION_LOCAL))
    parser.add_argument("--phenotype", default=str(TCGA_SKCM_PHENOTYPE_LOCAL))
    parser.add_argument("--survival", default=str(TCGA_SKCM_SURVIVAL_LOCAL))
    parser.add_argument(
        "--source",
        choices=["xena", "cbioportal"],
        default="xena",
        help="Bulk data source format. Use cbioportal for study files from the cBioPortal datahub.",
    )
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    REFS_DIR.mkdir(parents=True, exist_ok=True)
    rankings = pd.read_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t")
    top_tfs = rankings["tf"].head(args.top_n).tolist()

    expr_path = Path(args.expr)
    phenotype_path = Path(args.phenotype)
    survival_path = Path(args.survival)

    if args.source == "cbioportal":
        expr = read_cbioportal_table(expr_path)
        phenotype = read_cbioportal_table(phenotype_path)
        survival = read_cbioportal_table(survival_path)
        expr = expr.rename(columns={"Hugo_Symbol": "gene_symbol", "Entrez_Gene_Id": "entrez_gene_id"})
        tf_expr = expr.loc[expr["gene_symbol"].isin(top_tfs)].copy()
        sample_columns = [col for col in tf_expr.columns if col not in {"gene_symbol", "entrez_gene_id"}]
        melted = tf_expr.melt(
            id_vars=["gene_symbol", "entrez_gene_id"],
            value_vars=sample_columns,
            var_name="sample",
            value_name="log2_fpkm_uq_plus_1",
        )
        phenotype = phenotype.rename(columns={"SAMPLE_ID": "sample", "PATIENT_ID": "patient_id"})
        phenotype["disease_group"] = phenotype["SAMPLE_TYPE"].map(classify_cbioportal_sample_type)
    else:
        expr = read_table(expr_path if expr_path.exists() else TCGA_SKCM_EXPRESSION_URL)
        phenotype = read_table(phenotype_path if phenotype_path.exists() else TCGA_SKCM_PHENOTYPE_URL)
        survival = read_table(survival_path if survival_path.exists() else TCGA_SKCM_SURVIVAL_URL)

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

        phenotype_sample_col = phenotype.columns[0]
        disease_group_col = find_first_column(phenotype, TCGA_SKCM_PRIMARY_METASTATIC_CANDIDATES)
        phenotype = phenotype.rename(columns={phenotype_sample_col: "sample"})
        if disease_group_col is not None:
            phenotype["disease_group"] = phenotype[disease_group_col].map(classify_primary_vs_metastatic)
        else:
            phenotype["disease_group"] = None

    summary = (
        melted.groupby("gene_symbol", as_index=False)
        .agg(
            mean_bulk_expr=("log2_fpkm_uq_plus_1", "mean"),
            median_bulk_expr=("log2_fpkm_uq_plus_1", "median"),
            std_bulk_expr=("log2_fpkm_uq_plus_1", "std"),
            n_samples=("log2_fpkm_uq_plus_1", "count"),
        )
        .rename(columns={"gene_symbol": "tf"})
    )

    melted = melted.merge(phenotype[["sample", "disease_group"]], on="sample", how="left")
    grouped = (
        melted.dropna(subset=["disease_group"])
        .groupby(["gene_symbol", "disease_group"], as_index=False)
        .agg(
            mean_log2_fpkm_uq_plus_1=("log2_fpkm_uq_plus_1", "mean"),
            median_log2_fpkm_uq_plus_1=("log2_fpkm_uq_plus_1", "median"),
            n_samples=("log2_fpkm_uq_plus_1", "count"),
        )
        .rename(columns={"gene_symbol": "tf"})
    )

    phenotype.to_csv(RESULTS_DIR / "tcga_skcm_phenotype.tsv", sep="\t", index=False)
    survival.to_csv(RESULTS_DIR / "tcga_skcm_survival.tsv", sep="\t", index=False)
    melted.to_csv(RESULTS_DIR / "tcga_skcm_tf_expression_long.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS_DIR / "tcga_skcm_tf_expression_summary.tsv", sep="\t", index=False)
    grouped.to_csv(RESULTS_DIR / "tcga_skcm_tf_expression_by_group.tsv", sep="\t", index=False)

    plot_top_tfs(summary, FIGURES_DIR / "tcga_skcm_top_tf_expression.png")
    plot_primary_vs_metastatic(grouped, FIGURES_DIR / "tcga_skcm_primary_vs_metastatic.png")


if __name__ == "__main__":
    main()
