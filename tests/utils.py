import pandas as pd
import numpy as np


def parse_genome_file(filepath: str) -> pd.DataFrame:
    """Helper to parse a .genome file into a df"""
    df = pd.read_csv(filepath, sep=r"\s+", engine="python")

    assert df.columns.tolist() == [
        "FID1",
        "IID1",
        "FID2",
        "IID2",
        "RT",
        "EZ",
        "Z0",
        "Z1",
        "Z2",
        "PI_HAT",
        "PHE",
        "DST",
        "PPC",
        "RATIO",
    ]

    return df


def calculate_r_squared(truth, target):
    """Calculates the R^2 correlation between two lists of numbers."""
    if len(truth) != len(target):
        return 0.0

    correlation_matrix = np.corrcoef(np.array(truth), np.array(target))
    r = correlation_matrix[0, 1]

    # edge case, where truth or target is all 0s then correlation is not well defined=
    if np.isnan(r):
        r = 1.0 if all(x == 0 for x in truth) and all(x == 0 for x in target) else 0.0

    return r**2


def validate_genome_data(
    truth_filepath: str, target_filepath: str, r_squared_threshold: float = 0.99
):
    """Helper to validate that the target output matches the ground truth"""

    truth_data = parse_genome_file(truth_filepath)
    target_data = parse_genome_file(target_filepath)

    # 1. validate shapes match

    assert (
        truth_data.shape == target_data.shape
    ), f"Expected shape {truth_data.shape} but got {target_data.shape}"

    # 2. validate columns match

    assert (
        truth_data.columns.tolist() == target_data.columns.tolist()
    ), f"Expected columns {truth_data.columns.tolist()} but got {target_data.columns.tolist()}"

    # 3. validate [FID1, IID1, FID2, IID2] are in the same order in both dataframes

    assert (
        truth_data[["FID1", "IID1", "FID2", "IID2"]].values.tolist()
        == target_data[["FID1", "IID1", "FID2", "IID2"]].values.tolist()
    ), f"Expected FID1, IID1, FID2, IID2 order to match between truth and target, but they do not match"

    # 4. validate Z0, Z1, Z2, PI_HAT are similar enough (R^2 > r_squared_threshold)

    for column in ["Z0", "Z1", "Z2", "PI_HAT"]:
        r_squared = calculate_r_squared(truth_data[column], target_data[column])
        assert (
            r_squared > r_squared_threshold
        ), f"R^2 for column {column} is {r_squared}, which is below the threshold of {r_squared_threshold}"
