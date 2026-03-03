import logging
import os

import pandas as pd
from bed_reader import open_bed

import sys

logger = logging.getLogger("python_ibd")


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


def print_progress(
    name: str,
    sub_name: str = None,
    sub_current: int = None,
    sub_total: int = None,
    bar_width: int = 30,
    is_finished=False,
    print_all_iterations=False,
) -> None:
    """Write an inline progress bar to stdout using carriage return."""

    # only print every 0.1% of iterations by default to avoid overwhelming stdout, but print all iterations.
    if (
        not print_all_iterations
        and sub_current is not None
        and sub_total is not None
        and sub_total > 1000
        and sub_current % (sub_total // 1000) != 0
        and not is_finished
    ):
        return

    # Base output string
    output = f"\r{name}"

    if (sub_name or (sub_current and sub_total)) and not is_finished:
        pct = (sub_current + 1) / sub_total
        filled = int(bar_width * pct)
        bar = (
            "=" * filled
            + (">" if filled < bar_width else "")
            + " " * max(0, bar_width - filled - 1)
        )
        output += f" | {sub_name} [{bar}] {pct * 100:5.1f}%"

    # \033[K clears the rest of the line to prevent trailing artifact characters
    sys.stdout.write(f"{output}\033[K")
    sys.stdout.flush()

    if is_finished:
        sys.stdout.write("\n")
        sys.stdout.flush()


def write_genome_file(output_prefix: str, result: pd.DataFrame):
    """Helper to write a dataframe to a .genome file with the correct formatting"""

    # make file if does not exist, otherwise overwrite
    filepath = f"{output_prefix}.genome"

    dir_name = os.path.dirname(filepath)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    with open(filepath, "w") as f:
        f.write(result.to_string(justify="right", index=False))


def run_implementation(ibd_fn, input_prefix, out_prefix, implementation_name=None):
    setup_logger(out_prefix)
    implementation_name = implementation_name or "unknown"
    logger.debug(
        "IBD calculation recieved parameters input prefix: '%s' and output prefix: '%s'",
        input_prefix,
        out_prefix,
    )

    print_progress("Stage 1/5 Reading Input File...", is_finished=False)
    logger.debug("Started reading input file.")

    bed = open_bed(f"{input_prefix}.bed")
    individual_ids = bed.iid
    family_ids = bed.fid
    genotypes = bed.read()

    print_progress("Stage 1/5 Reading Input File... Done.", is_finished=True)
    logger.debug("Finished reading input file.")

    logger.debug(
        "Input files contain %d individuals, %d variants.",
        genotypes.shape[0],
        genotypes.shape[1],
    )

    ibd_results = ibd_fn(genotypes)

    result = pd.DataFrame(
        columns=[
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
        ],
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

    print_progress("Stage 5/5 Writing to output files...", is_finished=False)

    write_genome_file(out_prefix, result)

    print_progress(
        f"Stage 5/5 Writing to output files... Done.",
        is_finished=True,
    )

    logger.info("Output written to %s.genome and %s.log", out_prefix, out_prefix)
