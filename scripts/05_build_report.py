from __future__ import annotations

import pandas as pd

import _bootstrap  # noqa: F401
from config.config import RESULTS_DIR


def main() -> None:
    perturb = pd.read_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t")
    tcga = pd.read_csv(RESULTS_DIR / "tcga_skcm_tf_expression_summary.tsv", sep="\t")
    mode_values = set(perturb["mode"].dropna().unique())
    if mode_values == {"ranked_fallback"}:
        merged = perturb.merge(tcga, on="tf", how="left")
        merged["report_status"] = "diagnostic_only"
        merged["report_warning"] = (
            "Fallback perturbation values come from TF-IDF token deletion, not Geneformer. "
            "This table is for software inspection only and not for biological ranking."
        )
        merged = merged.sort_values(["mean_probability_drop", "tf"], ascending=[False, True])
        merged.to_csv(RESULTS_DIR / "candidate_tf_report.tsv", sep="\t", index=False)
        print(merged.head(20).to_string(index=False))
        return

    raise RuntimeError(
        "Unsupported perturbation mode for report generation. "
        "Add a mode-specific reporting path once official Geneformer outputs exist."
    )


if __name__ == "__main__":
    main()
