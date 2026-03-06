import sys
import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from tests.shared import parse_genome_file, calculate_r_squared

NUMERIC_COLS = ["Z0", "Z1", "Z2", "PI_HAT"]


def compare(truth_path: str, target_path: str):
    """
    Compare two .genome files and prints accuracy metrics (R², MAE, Max Error) for the Z0, Z1, Z2, and PI_HAT columns.
    """
    truth = parse_genome_file(truth_path)
    target = parse_genome_file(target_path)

    assert (
        truth.shape == target.shape
    ), f"Shape mismatch: truth {truth.shape} vs target {target.shape}"

    # Ensure the same ordering of pairs in both files
    id_cols = ["FID1", "IID1", "FID2", "IID2"]
    assert (
        truth[id_cols].values.tolist() == target[id_cols].values.tolist()
    ), "Individual pair ordering differs between files"

    n_pairs = len(truth)

    print(f"Comparing {n_pairs} individual pairs\n")

    rows = []
    for col in NUMERIC_COLS:
        t = np.array(truth[col], dtype=float)
        p = np.array(target[col], dtype=float)
        diff = np.abs(t - p)
        rows.append(
            {
                "Column": col,
                "R²": calculate_r_squared(t, p),
                "MAE": np.mean(diff),
                "Max Error": np.max(diff),
            }
        )

    all_truth = np.concatenate([np.array(truth[c], dtype=float) for c in NUMERIC_COLS])
    all_target = np.concatenate(
        [np.array(target[c], dtype=float) for c in NUMERIC_COLS]
    )
    all_diff = np.abs(all_truth - all_target)
    rows.append(
        {
            "Column": "Overall",
            "R²": calculate_r_squared(all_truth, all_target),
            "MAE": np.mean(all_diff),
            "Max Error": np.max(all_diff),
        }
    )

    df = pd.DataFrame(rows).set_index("Column")
    print(df.to_string())


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python3 benchmarking/accuracy.py <ground_truth.genome> <target.genome>",
            file=sys.stderr,
        )
        sys.exit(1)

    compare(sys.argv[1], sys.argv[2])
