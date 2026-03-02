import itertools
import logging
import math

import numpy as np
import pandas as pd

from src.utils import (
    print_progress,
    run_implementation,
)

logger = logging.getLogger("python_ibd")


def compute_ibd(genotypes: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the IBD estimates (Z0, Z1, Z2) for all individuals in the genotype matrix using the method of moments approach.
    This is the naive implementation that should match PLINK's output but is not optimized for performance.
    """

    NUM_VARIANTS = genotypes.shape[1]
    NUM_INDIVIDUALS = genotypes.shape[0]

    logger.info("Computing allele frequencies and global expected IBS counts...")
    logger.debug("  %d variants, %d individuals", NUM_VARIANTS, NUM_INDIVIDUALS)

    # Step 1. For each SNP we precompute allele frequencies:
    # - X <- Sum reference allele count across all individuals
    # - Y <- Sum alternate allele count across all samples
    # - T <- Total allele count across all samples (T = X + Y)
    # - p <- reference allele frequency (p = X / T)
    # - q <- alternate allele frequency (q = Y / T)
    # variant_stats stores this information: variant_stats[snp_index] = { "ref_count": ..., "alt_count": ... }

    variant_stats = compute_variant_stats(genotypes, NUM_VARIANTS)

    # Step 2. Compute average global expected counts of IBS states conditional on IBD states

    avg_e00, avg_e01, avg_e02, avg_e11, avg_e12 = compute_average_expected_counts(
        NUM_VARIANTS, variant_stats
    )

    # Step 3. compute per individual pair IBD estimation
    # Purcell et al. describes that the per pair N(I | Z) is computed only via variants where both individuals have non-missing genotypes:
    # - pair specific  N(I | Z) <- gloabl expected P(IBS | IBD) average * S, where S is the count of non-missing SNPs for the pair of individuals
    # Results stored in matrix `result` where columns are individual 1 index, individual 2 index, Z0, Z1, Z2

    pairs = list(itertools.combinations(range(NUM_INDIVIDUALS), 2))
    total_pairs = len(pairs)
    logger.info("Done. Computing pairwise IBD for %d pairs:", total_pairs)

    result = np.zeros((total_pairs, 5))
    for pair_inx, (individual_1_idx, individual_2_idx) in enumerate(pairs):
        print_progress(pair_inx, total_pairs)
        ibs_0_count = 0
        ibs_1_count = 0
        ibs_2_count = 0

        for snp_index in range(NUM_VARIANTS):
            g1 = genotypes[individual_1_idx, snp_index]
            g2 = genotypes[individual_2_idx, snp_index]

            if math.isnan(g1) or math.isnan(g2):
                continue

            diff = abs(g1 - g2)
            if diff == 2:
                ibs_0_count += 1  # opposite homozygotes: (0,2) or (2,0)
            elif diff == 0:
                ibs_2_count += 1  # same genotype: (0,0), (1,1) or (2,2)
            else:
                ibs_1_count += 1  # one allele shared: (0,1), (1,0), (1,2) or (2,1)

        S = ibs_0_count + ibs_1_count + ibs_2_count

        # We cannot compute IBD if individuals have no non-missing variants in common
        if S == 0:
            result[pair_inx, 0] = individual_1_idx
            result[pair_inx, 1] = individual_2_idx
            continue

        e00 = avg_e00 * S  # N(I=0 count | Z=0)
        e01 = avg_e01 * S  # N(I=1 count | Z=1)
        e02 = avg_e02 * S  # N(I=2 count | Z=2)
        e11 = avg_e11 * S  # N(I=1 count | Z=1)
        e12 = avg_e12 * S  # N(I=2 count | Z=1)
        e22 = 1.0 * S  # N(I=2 count | Z=2) = S

        z0 = ibs_0_count / e00 if e00 > 0 else 0.0
        z1 = (ibs_1_count - z0 * e01) / e11 if e11 > 0 else 0.0
        z2 = (ibs_2_count - z0 * e02 - z1 * e12) / e22 if e22 > 0 else 0.0

        # Step 4. Apply bounding procedure
        z0, z1, z2 = bind_z_values(z0, z1, z2)

        result[pair_inx, 0] = individual_1_idx
        result[pair_inx, 1] = individual_2_idx
        result[pair_inx, 2] = z0
        result[pair_inx, 3] = z1
        result[pair_inx, 4] = z2

    return result


def bind_z_values(z0, z1, z2):
    """Helper function to apply the bounding procedure described in Purcell et al. 2007 to ensure valid IBD estimates."""

    # Purcell et al. defines the following binding procedure to ensure valid probabilities:
    # - If any Z value exceeds 1, clamp it to 1 and set the other two to 0
    # - If any Z value is negative, set it to 0 and renormalize the other two to sum to 1

    # clamp to 1 and 0 if any exceed 1
    if z0 > 1:
        z0, z1, z2 = 1.0, 0.0, 0.0
    if z1 > 1:
        z1, z0, z2 = 1.0, 0.0, 0.0
    if z2 > 1:
        z2, z0, z1 = 1.0, 0.0, 0.0

    # renormalize if any are negative
    if z0 < 0.0:
        z0 = 0.0
        z1_z2_sum = z1 + z2
        if z1_z2_sum > 0:
            z1 /= z1_z2_sum
            z2 /= z1_z2_sum
    if z1 < 0.0:
        z1 = 0.0
        z0_z2_sum = z0 + z2
        if z0_z2_sum > 0:
            z0 /= z0_z2_sum
            z2 /= z0_z2_sum
    if z2 < 0.0:
        z2 = 0.0
        z0_z1_sum = z0 + z1
        if z0_z1_sum > 0:
            z0 /= z0_z1_sum
            z1 /= z0_z1_sum
    return z0, z1, z2


def compute_average_expected_counts(NUM_VARIANTS: int, variant_stats: dict) -> tuple:
    """Helper function to compute the global expected counts of IBS states conditional on IBD states, used in the method of moments IBD estimation."""

    # Purcell et al. describes the following:
    # - N(I = i | Z = z) <- The global expected count of SNPs with IBS state I = i conditional on IBD state Z = z.
    # - N(I = i | Z = z) <- Computed via ∑_m P(IBS = i | IBD = z) over all SNPs m
    # - P(IBS = i | IBD = z) can be computed via Table 1. of Purcell et al. 2007
    # Plink 1.9 computes these global expected via the following notation:
    # - sum_e00 <- N(I=0 | Z=0)
    # - sum_e01 <- N(I=1 | Z=0)
    # - sum_e02 <- N(I=2 | Z=0)
    # - sum_e11 <- N(I=1 | Z=1)
    # - sum_e12 <- N(I=2 | Z=1)

    sum_e00 = 0.0
    sum_e01 = 0.0
    sum_e02 = 0.0
    sum_e11 = 0.0
    sum_e12 = 0.0
    cnt_poly = 0  # count of polymorphic SNPs, used by PLINK 1.9 (ibd_prect) in per-indivual pair caluclations

    for variant_id in range(NUM_VARIANTS):
        # attempt to match Purcell et al. then PLINK 1.9's naming conventions as much as possible, see https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4849-L4866
        X = variant_stats[variant_id]["ref_count"]  # ref count
        Y = variant_stats[variant_id]["alt_count"]  # alt count
        T = X + Y  # total variant count

        # Skip monomorphic SNPs (X == 0 or Y == 0) and SNPs with low count (T < 4), see https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4846C2-L4846C56
        if X == 0 or Y == 0 or T < 4:
            continue

        p = X / T  # ref allele frequency
        q = Y / T  # alt allele frequency
        dpp_sq = p**2
        dqq_sq = q**2
        dxx1 = (X - 1) / X if X > 0 else 0.0
        dxx2 = (X - 1) * (X - 2) / (X**2) if X > 1 else 0.0
        dyy1 = (Y - 1) / Y if Y > 0 else 0.0
        dyy2 = (Y - 1) * (Y - 2) / (Y**2) if Y > 1 else 0.0
        num_allelesf2 = T * T / ((T - 1) * (T - 2))
        num_allelesf3 = num_allelesf2 * T / (T - 3)

        sum_e00 += 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3
        sum_e01 += 4 * p * q * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2)
        sum_e02 += num_allelesf3 * (
            dqq_sq * dqq_sq * dyy2 * (Y - 3) / Y
            + dpp_sq * dpp_sq * dxx2 * (X - 3) / X
            + 4 * dpp_sq * dqq_sq * dxx1 * dyy1
        )
        sum_e11 += 2 * p * q * num_allelesf2 * (p * dxx1 + q * dyy1)
        sum_e12 += num_allelesf2 * (
            dpp_sq * p * dxx2
            + dqq_sq * q * dyy2
            + dpp_sq * q * dxx1
            + p * dqq_sq * dyy1
        )
        cnt_poly += 1

    # Compute average expected counts of IBS states, see https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4888

    avg_e00 = sum_e00 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e01 = sum_e01 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e02 = sum_e02 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e11 = sum_e11 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e12 = sum_e12 / cnt_poly if cnt_poly > 0 else 0.0
    return avg_e00, avg_e01, avg_e02, avg_e11, avg_e12


def compute_variant_stats(genotypes: np.ndarray, NUM_VARIANTS: int) -> dict:
    """Helper function to compute allele counts for each variant, used in the method of moments IBD estimation."""
    variant_stats = {}
    for variant_id in range(NUM_VARIANTS):
        variant_stats[variant_id] = {"ref_count": 0, "alt_count": 0}
        for sample_id in range(genotypes.shape[0]):
            genotype_for_sample = genotypes[sample_id, variant_id]
            if not math.isnan(genotype_for_sample):
                if genotype_for_sample == 0:
                    variant_stats[variant_id]["ref_count"] += 2
                elif genotype_for_sample == 1:
                    variant_stats[variant_id]["ref_count"] += 1
                    variant_stats[variant_id]["alt_count"] += 1
                elif genotype_for_sample == 2:
                    variant_stats[variant_id]["alt_count"] += 2
    return variant_stats


def run_naive(input_prefix, out_prefix):
    run_implementation(
        compute_ibd, input_prefix, out_prefix, implementation_name="naive"
    )
