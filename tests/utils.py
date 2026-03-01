import pandas as pd
import numpy as np
import tempfile
import os


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
    """R²^2 as coefficient of determination: 1 - SS_res / SS_tot."""
    truth = np.array(truth, dtype=float)
    target = np.array(target, dtype=float)

    if len(truth) != len(target) or len(truth) == 0:
        return 0.0

    ss_res = np.sum((truth - target) ** 2)
    ss_tot = np.sum((truth - np.mean(truth)) ** 2)

    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0

    return 1.0 - ss_res / ss_tot


def assert_genome_data(
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


def get_tmp_output_prefix(prefix: str = "test_output") -> str:
    """Helper to generate a temporary output prefix for testing"""
    tmpdir = tempfile.mkdtemp()
    return f"{tmpdir}/{prefix}"


def assert_log_file_exists(out_prefix):
    assert os.path.exists(
        f"{out_prefix}.log"
    ), "Expected log file to be generated but it does not exist"
