from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import _bootstrap  # noqa: F401
from config.config import RESULTS_DIR, SEED_TFS


def score_embedding_shift(oracle) -> float:
    if hasattr(oracle, "delta_embedding"):
        delta = np.asarray(oracle.delta_embedding)
        return float(np.linalg.norm(delta, axis=1).mean())
    if hasattr(oracle, "delta_embedding_random"):
        delta = np.asarray(oracle.delta_embedding_random)
        return float(np.linalg.norm(delta, axis=1).mean())
    return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tf",
        action="append",
        dest="tfs",
        default=[],
        help="TF to perturb. Can be supplied multiple times.",
    )
    parser.add_argument("--oracle", default=str(RESULTS_DIR / "celloracle_oracle.hdf5"))
    args = parser.parse_args()

    try:
        import celloracle as co
    except ImportError as exc:
        raise RuntimeError(
            "CellOracle is not available in this environment. "
            "Use a Linux or Docker-backed environment that can install CellOracle cleanly."
        ) from exc

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    oracle = co.load_hdf5(args.oracle)
    tfs = args.tfs or SEED_TFS

    rows = []
    for tf in tfs:
        oracle_sim = oracle.copy()
        oracle_sim.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        oracle_sim.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1.0)
        oracle_sim.calculate_embedding_shift(sigma_corr=0.05)
        rows.append(
            {
                "tf": tf,
                "celloracle_embedding_shift_norm": score_embedding_shift(oracle_sim),
                "mode": "celloracle",
            }
        )

    out = pd.DataFrame(rows).sort_values("celloracle_embedding_shift_norm", ascending=False)
    out.to_csv(RESULTS_DIR / "celloracle_tf_perturbation_ranking.tsv", sep="\t", index=False)
    with open(RESULTS_DIR / "celloracle_run_metadata.json", "w") as handle:
        json.dump({"mode": "celloracle", "tested_tfs": tfs}, handle, indent=2)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
