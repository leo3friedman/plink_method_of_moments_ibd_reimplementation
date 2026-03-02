import logging
import os

import pandas as pd
from bed_reader import open_bed

logger = logging.getLogger("python_ibd")

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


def setup_logger(output_prefix: str) -> logging.Logger:
    """Configure the 'python_ibd' logger to write to {output_prefix}.log and stdout."""
    filepath = f"{output_prefix}.log"
    dir_name = os.path.dirname(filepath)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    logger = logging.getLogger("python_ibd")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False  # prevent double-printing via root logger

    file_fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    fh = logging.FileHandler(filepath, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)
    logger.addHandler(fh)

    # stdout only shows INFO+, with no timestamp — just the message
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    return logger


def print_progress(current: int, total: int, bar_width: int = 40) -> None:
    """Write an inline progress bar to stdout using carriage return."""
    import sys

    pct = (current + 1) / total
    filled = int(bar_width * pct)
    bar = (
        "=" * filled
        + (">" if filled < bar_width else "")
        + " " * max(0, bar_width - filled - 1)
    )
    sys.stdout.write(f"\r  [{bar}] {pct * 100:5.1f}%  ({current + 1}/{total})")
    sys.stdout.flush()
    if current == total - 1:
        sys.stdout.write("\n")
        sys.stdout.flush()


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


def run_implementation(ibd_fn, input_prefix, out_prefix, implementation_name=None):
    setup_logger(out_prefix)
    implementation_name = implementation_name or "unknown"
    logger.info(
        "Running %s IBD | input=%s | out=%s",
        implementation_name,
        input_prefix,
        out_prefix,
    )

    logger.info("Reading %s.bed...", input_prefix)
    bed = open_bed(f"{input_prefix}.bed")
    individual_ids = bed.iid
    family_ids = bed.fid
    genotypes = bed.read()
    logger.info(
        "Done. %d individuals, %d variants.",
        genotypes.shape[0],
        genotypes.shape[1],
    )

    ibd_results = ibd_fn(genotypes)

    logger.debug("Building result dataframe and writing output...")
    result = pd.DataFrame(
        columns=GENOME_COLUMNS,
        data=[
            {
                "FID1": family_ids[int(sample_idx1)],
                "IID1": individual_ids[int(sample_idx1)],
                "FID2": family_ids[int(sample_idx2)],
                "IID2": individual_ids[int(sample_idx2)],
                "RT": "UN",
                "EZ": "NA",
                "Z0": Z0,
                "Z1": Z1,
                "Z2": Z2,
                "PI_HAT": Z1 / 2 + Z2,
                "PHE": -1,
                "DST": 0,
                "PPC": 0,
                "RATIO": 0,
            }
            for (sample_idx1, sample_idx2, Z0, Z1, Z2) in list(ibd_results)
        ],
    )

    write_genome_file(out_prefix, result)
    logger.info("Done. Output written to %s.genome and %s.log", out_prefix, out_prefix)
