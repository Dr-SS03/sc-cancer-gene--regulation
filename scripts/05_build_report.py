from __future__ import annotations

import pandas as pd

import _bootstrap  # noqa: F401
from config.config import RESULTS_DIR


def main() -> None:
    perturb = pd.read_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t")
    tcga = pd.read_csv(RESULTS_DIR / "tcga_skcm_tf_expression_summary.tsv", sep="\t")
    merged = perturb.merge(tcga, on="tf", how="left")
    merged["priority_score"] = merged["mean_probability_drop"] * merged["mean_bulk_expr"].fillna(0.0)
    merged = merged.sort_values("priority_score", ascending=False)
    merged.to_csv(RESULTS_DIR / "candidate_tf_report.tsv", sep="\t", index=False)
    print(merged.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
