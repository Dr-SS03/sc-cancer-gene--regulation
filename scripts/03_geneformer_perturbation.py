from __future__ import annotations

import argparse
import json

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfVectorizer

import _bootstrap  # noqa: F401
from config.config import RESULTS_DIR, SEED_TFS


def compute_geneformer_embedding_matrix(
    ranked_cells: pd.DataFrame, vectorizer: TfidfVectorizer | None = None
) -> tuple[np.ndarray, TfidfVectorizer]:
    try:
        from geneformer import EmbExtractor
    except ImportError as exc:
        raise RuntimeError(
            "Geneformer is not installed. Create the conda environment from environment.yml first."
        ) from exc

    raise NotImplementedError(
        "Official Geneformer mode requires tokenizer-compliant raw-count inputs. "
        "Per Geneformer docs, supply raw-count .h5ad/.loom with ensembl_id and n_counts, "
        "then wire TranscriptomeTokenizer -> EmbExtractor/InSilicoPerturber here."
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


def compute_ranked_fallback_embeddings(
    ranked_cells: pd.DataFrame, vectorizer: TfidfVectorizer | None = None
) -> tuple[np.ndarray, TfidfVectorizer]:
    corpus = ranked_cells["ranked_genes"].fillna("")
    if vectorizer is None:
        vectorizer = TfidfVectorizer(token_pattern=r"(?u)\b[A-Za-z0-9\-\.]+\b", max_features=4096)
        matrix = vectorizer.fit_transform(corpus)
    else:
        matrix = vectorizer.transform(corpus)
    return matrix.toarray(), vectorizer


def run_ranked_fallback(ranked: pd.DataFrame) -> dict[str, float | list[str]]:
    malignant_mask = ranked["malignant"].astype(str).eq("2").to_numpy()
    embeddings, vectorizer = compute_ranked_fallback_embeddings(ranked)
    clf, auc = train_malignancy_classifier(embeddings, malignant_mask.astype(int))
    joblib.dump(clf, RESULTS_DIR / "malignancy_classifier.joblib")
    joblib.dump(vectorizer, RESULTS_DIR / "ranked_fallback_vectorizer.joblib")

    perturbation_rows = []
    malignant_cells = ranked.loc[malignant_mask].copy()
    base_embeddings, _ = compute_ranked_fallback_embeddings(malignant_cells, vectorizer)
    base_scores = clf.predict_proba(base_embeddings)[:, 1]

    for tf in SEED_TFS:
        perturbed = malignant_cells.copy()
        perturbed["ranked_genes"] = perturbed["ranked_genes"].map(lambda value: perturb_ranked_genes(value, tf))
        perturbed_embeddings, _ = compute_ranked_fallback_embeddings(perturbed, vectorizer)
        perturbed_scores = clf.predict_proba(perturbed_embeddings)[:, 1]
        delta = base_scores - perturbed_scores
        perturbation_rows.append(
            {
                "tf": tf,
                "n_cells_scored": int(len(delta)),
                "mean_probability_drop": float(np.mean(delta)),
                "median_probability_drop": float(np.median(delta)),
                "mode": "ranked_fallback",
                "interpretation": "diagnostic_only",
            }
        )

    out = pd.DataFrame(perturbation_rows).sort_values("mean_probability_drop", ascending=False)
    out.to_csv(RESULTS_DIR / "geneformer_tf_perturbation_ranking.tsv", sep="\t", index=False)
    return {
        "classifier_auc": auc,
        "mode": "ranked_fallback",
        "seed_tfs": SEED_TFS,
        "status": "diagnostic_only",
        "warning": (
            "ranked_fallback is a TF-IDF baseline, not Geneformer. "
            "Do not interpret the perturbation magnitudes biologically."
        ),
    }


def run_official_geneformer(ranked: pd.DataFrame) -> dict[str, float | list[str]]:
    malignant_mask = ranked["malignant"].astype(str).eq("2").to_numpy()
    embeddings, _ = compute_geneformer_embedding_matrix(ranked)
    clf, auc = train_malignancy_classifier(embeddings, malignant_mask.astype(int))
    joblib.dump(clf, RESULTS_DIR / "malignancy_classifier.joblib")
    return {"classifier_auc": auc, "mode": "official", "seed_tfs": SEED_TFS}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["official", "ranked_fallback"],
        default="ranked_fallback",
        help="Choose official Geneformer mode for tokenizer-compliant raw-count inputs or ranked_fallback for the staged Tirosh matrix.",
    )
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ranked = pd.read_parquet(RESULTS_DIR / "geneformer_ranked_cells.parquet")

    with open(RESULTS_DIR / "geneformer_run_metadata.json", "w") as handle:
        if args.mode == "official":
            json.dump(run_official_geneformer(ranked), handle, indent=2)
        else:
            json.dump(run_ranked_fallback(ranked), handle, indent=2)


if __name__ == "__main__":
    main()
