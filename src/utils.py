import pandas as pd
import os

GENOME_COLUMNS = [
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


def write_genome_file(output_prefix: str, result: pd.DataFrame):
    """Helper to write a dataframe to a .genome file with the correct formatting"""

    # make file if does not exist, otherwise overwrite
    filepath = f"{output_prefix}.genome"
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, "w") as f:
        f.write(result.to_string(justify="right", index=False))


def write_log_file(output_prefix: str, content: str):
    """Helper to write content to a log file"""

    # make file if does not exist, otherwise overwrite
    filepath = f"{output_prefix}.log"
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, "w") as f:
        f.write(content)
