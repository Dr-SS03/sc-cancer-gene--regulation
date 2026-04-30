from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import urllib.request
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401
from config.config import GSE72056_URL, PROCESSED_DIR, RAW_DIR


def ensure_raw_dataset(local_gz: str | None) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    target = RAW_DIR / "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
    if local_gz:
        src = Path(local_gz).expanduser().resolve()
        if src != target:
            shutil.copy2(src, target)
    elif not target.exists():
        print(f"Downloading {GSE72056_URL} -> {target}")
        urllib.request.urlretrieve(GSE72056_URL, target)
    return target


def split_tirosh_matrix(gz_path: Path) -> dict[str, int]:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    txt_path = PROCESSED_DIR / "GSE72056_melanoma_single_cell_revised_v2.txt"
    with gzip.open(gz_path, "rb") as src, open(txt_path, "wb") as dst:
        shutil.copyfileobj(src, dst)

    metadata_path = PROCESSED_DIR / "cell_metadata.tsv"
    expression_path = PROCESSED_DIR / "expression_matrix.tsv"

    with open(txt_path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        tumor = next(reader)
        malignant = next(reader)
        celltype = next(reader)

        with open(metadata_path, "w", newline="") as meta_out:
            writer = csv.writer(meta_out, delimiter="\t")
            writer.writerow(["cell_id", "tumor", "malignant", "non_malignant_cell_type"])
            for idx, cell_id in enumerate(header[1:]):
                writer.writerow(
                    [cell_id, tumor[idx + 1], malignant[idx + 1], celltype[idx + 1]]
                )

        with open(expression_path, "w", newline="") as expr_out:
            writer = csv.writer(expr_out, delimiter="\t")
            writer.writerow(header)
            for row in reader:
                writer.writerow(row)

    malignant_counter = Counter(malignant[1:])
    summary = {
        "n_cells": len(header) - 1,
        "n_tumors": len(set(tumor[1:])),
        "n_malignant": malignant_counter.get("2", 0),
        "n_non_malignant": malignant_counter.get("1", 0),
        "n_unresolved": malignant_counter.get("0", 0),
    }
    with open(PROCESSED_DIR / "dataset_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--local-gz",
        default=None,
        help="Optional path to an already-downloaded GSE72056 gzip file.",
    )
    args = parser.parse_args()

    gz_path = ensure_raw_dataset(args.local_gz)
    summary = split_tirosh_matrix(gz_path)
    print(json.dumps({"staged_gz": str(gz_path), **summary}, indent=2))


if __name__ == "__main__":
    main()
