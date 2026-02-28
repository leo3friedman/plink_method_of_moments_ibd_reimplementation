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


def write_genome_file(df: pd.DataFrame, filepath: str):
    """Helper to write a dataframe to a .genome file with the correct formatting"""

    # make file if does not exist, otherwise overwrite
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, "w") as f:
        f.write(df.to_string(justify="right", index=False))
