import pandas as pd
import numpy as np
import math

from src.naive import (
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

    logger.debug(
        "Computing IBD estimates for %d variants and %d individuals...",
        NUM_VARIANTS,
        NUM_INDIVIDUALS,
    )

    logger.info("Computing allele frequencies...")

    # Step 1. For each SNP we precompute allele frequencies. See naive.py for more implementation notes.

    variant_stats = compute_variant_stats(genotypes)

    logger.debug("Done. Finished computing allele frequencies.")

    # Step 2. Compute average global expected counts of IBS states conditional on IBD states. See naive.py for moreimplementation notes.

    logger.info("Computing global expected IBS counts...")

    avg_e00, avg_e01, avg_e02, avg_e11, avg_e12 = compute_average_expected_counts(
        variant_stats
    )

    logger.debug("Done. Finished computing global expected IBS counts.")

    # Step 3. compute per individual pair IBD estimation. See naive.py for more implementation notes.

    total_pairs = math.comb(NUM_INDIVIDUALS, 2)
    logger.info("Computing pairwise IBD for %d pairs...", total_pairs)

    # is0[i,m] = 1 if individual i has genotype 0 at SNP m. Shape: (n_individuals × n_variants)
    # NaN genotypes evaluate to False for all is0, is1, and is2, so they are effectively ignored in the downstream IBS counts.
    is0 = (genotypes == 0).astype(np.float32)
    is1 = (genotypes == 1).astype(np.float32)
    is2 = (genotypes == 2).astype(np.float32)

    # valid[i, m] = 1 if individual i has a non-missing genotype at SNP m. Shape: (n_individuals × n.variants)
    valid = (~np.isnan(genotypes)).astype(np.float32)

    # ibs_2[i, j] = count of IBS 2 SNPs between individuals i and j. Shape: (n_individuals × n_individuals)
    # Note: (is0 @ is0.T)[i, j] is the count of SNPs where individuals i and j both have genotype 0. Similarily for is1 and is2, so summing gives the total count of matching genotypes.
    ibs_2 = is0 @ is0.T + is1 @ is1.T + is2 @ is2.T

    # ibs_0[i, j] = count of IBS 0 SNPs between individuals i and j. Shape: (n_individuals × n_individuals)
    # Note (is0 @ is2.T)[i, j] is the count of SNPs where individual i has genotype 0 and individual j has genotype 2. Similarily for (is2 @ is0.T), so summing gives the total count of IBS 0 SNPs.
    ibs_0 = is0 @ is2.T + is2 @ is0.T

    # S[i,j] = number of SNPs where both individuals have non-missing genotypes. Shape: (n_individuals × n_individuals)
    S = valid @ valid.T

    # ibs_1[i,j] = count of IBS 1 SNPs between individuals i and j. Shape: (n_individuals × n_individuals)
    # Can be computed as all non missing snps that is are not ibs_2 or ibs_0
    ibs_1 = S - ibs_2 - ibs_0

    # The indices of the upper triangle of the above (n_individuals × n_individuals) matrices
    # Each of these matrices are symmetric, so we don't distiguish between pair (i, j) and pair (j, i).
    # The diagonal of these matrices corresponds to self comparisons (i, i) which we also ignore.
    # (ind_i[i], ind_j[i]) gives the i'th pair of individuals for which we will compute IBD estimates.
    idx_i, idx_j = np.triu_indices(NUM_INDIVIDUALS, k=1)

    # Restrict the ibs matrices to only the upper right triangle pairs.
    # Shape: (n_pairs, 1), where ibs0[i] is the count of IBS 0 SNPs for pair i (i.e. pair (ind_i[i], ind_j[i]))
    ibs0 = ibs_0[idx_i, idx_j]
    ibs1 = ibs_1[idx_i, idx_j]
    ibs2 = ibs_2[idx_i, idx_j]
    S = S[idx_i, idx_j]

    # Compute expected IBS counts, see naive.py for more implementation notes.
    # Shape: (n_pairs, 1)
    e00 = avg_e00 * S  # N(I=0 count | Z=0)
    e01 = avg_e01 * S  # N(I=1 count | Z=1)
    e02 = avg_e02 * S  # N(I=2 count | Z=2)
    e11 = avg_e11 * S  # N(I=1 count | Z=1)
    e12 = avg_e12 * S  # N(I=2 count | Z=1)
    e22 = 1.0 * S  # N(I=2 count | Z=2) = S

    # Compute IBD estimates using method of moments equations, see naive.py for more implementation notes.
    # Shape: (n_pairs, 1)
    z0 = np.where(e00 > 0, ibs0 / e00, 0.0)
    z1 = np.where(e11 > 0, (ibs1 - z0 * e01) / e11, 0.0)
    z2 = np.where(e22 > 0, (ibs2 - z0 * e02 - z1 * e12) / e22, 0.0)

    logger.debug("Done. Finished computing pairwise IBD estimates.")

    # Step 4. Apply bounding procedure, see naive.py for more implementation notes.

    logger.info("Applying bounding procedure to IBD estimates...")

    for pair_idx in range(total_pairs):
        z0[pair_idx], z1[pair_idx], z2[pair_idx] = bind_z_values(
            z0[pair_idx], z1[pair_idx], z2[pair_idx]
        )

    logger.debug("Done. Finished applying bounding procedure.")

    result = np.zeros((total_pairs, 5))
    result[:, 0] = idx_i  # individual 1 index
    result[:, 1] = idx_j  # individual 2 index
    result[:, 2] = z0
    result[:, 3] = z1
    result[:, 4] = z2

    return result


def compute_variant_stats(genotypes: np.ndarray):
    """Helper function to compute allele frequencies and other variant stats for each SNP"""

    non_missing = ~np.isnan(genotypes)
    alt_counts = np.nansum(genotypes, axis=0).astype(int)  # sum down columns
    ref_counts = (2 * non_missing.sum(axis=0) - alt_counts).astype(int)
    return {
        v: {"ref_count": int(ref_counts[v]), "alt_count": int(alt_counts[v])}
        for v in range(genotypes.shape[1])
    }


def run_optimized(input_prefix, out_prefix):
    run_implementation(
        compute_ibd, input_prefix, out_prefix, implementation_name="optimized"
    )
