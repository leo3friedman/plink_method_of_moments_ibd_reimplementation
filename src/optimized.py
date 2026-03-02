import pandas as pd
import numpy as np

import math
import itertools
from src.naive import (
    compute_variant_stats,
    compute_average_expected_counts,
    bind_z_values,
)

import logging
from src.utils import print_progress, run_implementation

logger = logging.getLogger("python_ibd")


def compute_ibd(genotypes: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the IBD estimates (Z0, Z1, Z2) for all individuals in the genotype matrix using the method of moments approach.
    This is the optimized implementation which should match PLINK's output and should outperform the implementation in naive.py on large datasets.
    """

    NUM_VARIANTS = genotypes.shape[1]
    NUM_INDIVIDUALS = genotypes.shape[0]

    logger.info("Computing allele frequencies and global expected IBS counts...")
    logger.debug("  %d variants, %d individuals", NUM_VARIANTS, NUM_INDIVIDUALS)

    # Step 1. For each SNP we precompute allele frequencies. See naive.py for implementation notes.

    variant_stats = compute_variant_stats(genotypes, NUM_VARIANTS)

    # Step 2. Compute average global expected counts of IBS states conditional on IBD states. See naive.py for implementation notes.

    avg_e00, avg_e01, avg_e02, avg_e11, avg_e12 = compute_average_expected_counts(
        NUM_VARIANTS, variant_stats
    )

    # TODO: Step 3. compute per individual pair IBD estimation

    pairs = list(itertools.combinations(range(NUM_INDIVIDUALS), 2))
    total_pairs = len(pairs)
    logger.info("Computing pairwise IBD for %d pairs:", total_pairs)
    result = np.zeros((total_pairs, 5))
    result[:, 0] = [individual_1_idx for (individual_1_idx, _) in pairs]
    result[:, 1] = [individual_2_idx for (_, individual_2_idx) in pairs]

    return result


def run_optimized(input_prefix, out_prefix):
    run_implementation(
        compute_ibd, input_prefix, out_prefix, implementation_name="optimized"
    )
