from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

import _bootstrap  # noqa: F401
from config.config import RESULTS_DIR, SEED_TFS


def compute_geneformer_embedding_matrix(ranked_cells: pd.DataFrame) -> np.ndarray:
    try:
        from geneformer import EmbExtractor
    except ImportError as exc:
        raise RuntimeError(
            "Geneformer is not installed. Create the conda environment from environment.yml first."
        ) from exc

    # This script deliberately keeps the Geneformer call isolated in one place so the
    # repo can be inspected without importing the heavy model stack everywhere else.
    extractor = EmbExtractor(
        model_type="Pretrained",
        emb_mode="cell",
        summary_stat="mean",
    )
    # Geneformer expects tokenized inputs on disk in its own format. We keep the current
    # repo output simple and explicit, so adapt this block if you choose to use the
    # repository's tokenizer utilities directly.
    raise NotImplementedError(
        "Adapt EmbExtractor / tokenizer calls here to your installed Geneformer version. "
        "The surrounding classifier and ranking logic are ready once embeddings are produced."
    )


def train_malignancy_classifier(embeddings: np.ndarray, labels: np.ndarray) -> tuple[LogisticRegression, float]:
    x_train, x_test, y_train, y_test = train_test_split(
        embeddings,
        labels,
        test_size=0.2,
        random_state=7,
        stratify=labels,
    )
    model = LogisticRegression(max_iter=5000)
    model.fit(x_train, y_train)
    auc = roc_auc_score(y_test, model.predict_proba(x_test)[:, 1])
    return model, float(auc)


def perturb_ranked_genes(ranked_genes: str, tf: str) -> str:
    genes = ranked_genes.split()
    genes = [gene for gene in genes if gene != tf]
    return " ".join(genes)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ranked = pd.read_parquet(RESULTS_DIR / "geneformer_ranked_cells.parquet")
    malignant_mask = ranked["malignant"].astype(str).eq("2").to_numpy()

    embeddings = compute_geneformer_embedding_matrix(ranked)
    clf, auc = train_malignancy_classifier(embeddings, malignant_mask.astype(int))
    joblib.dump(clf, RESULTS_DIR / "malignancy_classifier.joblib")

    perturbation_rows = []
    malignant_cells = ranked.loc[malignant_mask].copy()
    base_scores = clf.predict_proba(embeddings[malignant_mask])[:, 1]

    for tf in SEED_TFS:
        perturbed = malignant_cells.copy()
        perturbed["ranked_genes"] = perturbed["ranked_genes"].map(lambda value: perturb_ranked_genes(value, tf))
        try:
            perturbed_embeddings = compute_geneformer_embedding_matrix(perturbed)
        except NotImplementedError:
            break
        perturbed_scores = clf.predict_proba(perturbed_embeddings)[:, 1]
        delta = base_scores - perturbed_scores
        perturbation_rows.append(
            {
                "tf": tf,
                "n_cells_scored": int(len(delta)),
                "mean_probability_drop": float(np.mean(delta)),
                "median_probability_drop": float(np.median(delta)),
            }
        )

    pd.DataFrame(perturbation_rows).sort_values(
        "mean_probability_drop", ascending=False
    ).to_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t", index=False)

    with open(RESULTS_DIR / "geneformer_run_metadata.json", "w") as handle:
        json.dump({"classifier_auc": auc, "seed_tfs": SEED_TFS}, handle, indent=2)


if __name__ == "__main__":
    main()
