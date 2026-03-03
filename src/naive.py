import itertools
import logging
import math

import numpy as np
import pandas as pd

from src.utils import (
    compute_expected_ibs,
    print_progress,
    run_implementation,
)

logger = logging.getLogger("python_ibd")


def compute_ibd(genotypes: pd.DataFrame) -> pd.DataFrame:
    """Computes the IBD estimates (Z0, Z1, Z2) for all individuals in the genotype matrix using the method of moments approach.

    This is the naive (loop-based) implementation which is not optimized for performance.
    """

    NUM_VARIANTS = genotypes.shape[1]
    NUM_INDIVIDUALS = genotypes.shape[0]

    # Step 1. Compute reference and alternate allele counts for each SNP.

    variant_stats = (
        {}
    )  # variant_stats[snp_index] = { "ref_count": ..., "alt_count": ... }

    for variant_id in range(genotypes.shape[1]):
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

        is_finished = variant_id == genotypes.shape[1] - 1
        print_progress(
            f"Stage 2/5 Computing Allele Frequencies... {'Done.' if is_finished else ''}",
            sub_name="Progress",
            sub_current=variant_id,
            sub_total=genotypes.shape[1],
            is_finished=is_finished,
        )

    # Step 2. Compute average global expected counts of IBS states conditional on IBD states.

    logger.debug("Started computing global expected IBS counts...")
    sum_e00 = 0.0  # N(I=0 | Z=0) <- global expected count of SNPs with IBS state 0 given IBD state 0
    sum_e01 = 0.0  # N(I=1 | Z=0)
    sum_e02 = 0.0  # N(I=2 | Z=0)
    sum_e11 = 0.0  # N(I=1 | Z=1)
    sum_e12 = 0.0  # N(I=2 | Z=1)
    cnt_poly = 0  # count of polymorphic SNPs with T > 3

    for i, variant_stat in enumerate(variant_stats.values()):
        is_finished = i == len(variant_stats) - 1
        print_progress(
            f"Stage 3/5 Computing IBS Estimates... {'Done.' if is_finished else ''}",
            sub_name="Progress",
            sub_current=i,
            sub_total=len(variant_stats),
            is_finished=is_finished,
        )

        X = variant_stat["ref_count"]
        Y = variant_stat["alt_count"]
        T = X + Y

        # Skip monomorphic SNPs (X == 0 or Y == 0) and SNPs with low count (T < 4)
        if X == 0 or Y == 0 or T < 4:
            continue

        e00, e01, e02, e11, e12 = compute_expected_ibs(X, Y)
        sum_e00 += e00
        sum_e01 += e01
        sum_e02 += e02
        sum_e11 += e11
        sum_e12 += e12
        cnt_poly += 1

    avg_e00 = sum_e00 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e01 = sum_e01 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e02 = sum_e02 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e11 = sum_e11 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e12 = sum_e12 / cnt_poly if cnt_poly > 0 else 0.0
    logger.debug("Finished computing global expected IBS counts.")

    # Step 3. Compute per individual pair IBD estimation.

    pairs = list(itertools.combinations(range(NUM_INDIVIDUALS), 2))
    total_pairs = len(pairs)

    result = np.zeros(
        (total_pairs, 5)
    )  # columns: [individual_1_idx, individual_2_idx, Z0, Z1, Z2]

    for pair_inx, (individual_1_idx, individual_2_idx) in enumerate(pairs):
        ibs_0_count = 0
        ibs_1_count = 0
        ibs_2_count = 0

        for snp_index in range(NUM_VARIANTS):
            g1 = genotypes[individual_1_idx, snp_index]
            g2 = genotypes[individual_2_idx, snp_index]

            # We only consider SNPs were both individuals have non-missing genotypes
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

        # We cannot compute IBD if the two individuals do not share any SNPs with non-missing genotypes
        if S == 0:
            result[pair_inx, 0] = individual_1_idx
            result[pair_inx, 1] = individual_2_idx
            continue

        # Per-pair expected IBS counts = (average global expected IBS count) * (number of shared variants S between the pair)
        e00 = avg_e00 * S  # N(I=0 | Z=0)
        e01 = avg_e01 * S  # N(I=0 | Z=1)
        e02 = avg_e02 * S  # N(I=2 | Z=0)
        e11 = avg_e11 * S  # N(I=1 | Z=1)
        e12 = avg_e12 * S  # N(I=2 | Z=1)
        e22 = 1.0 * S  # N(I=2 | Z=2) = S

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

        is_finished = pair_inx == total_pairs - 1
        print_progress(
            f"Stage 4/5 Computing pairwise IBD... {'Done.' if is_finished else ''}",
            sub_name=f"Progress",
            sub_current=pair_inx,
            sub_total=total_pairs,
            is_finished=is_finished,
        )

    return result


def bind_z_values(z0, z1, z2):
    """Implements the bounding procedure described in Purcell et al. 2007 to ensure valid IBD estimates."""

    # If any Z value exceeds 1, clamp it to 1 and set the other two to 0
    if z0 > 1:
        z0, z1, z2 = 1.0, 0.0, 0.0
    if z1 > 1:
        z1, z0, z2 = 1.0, 0.0, 0.0
    if z2 > 1:
        z2, z0, z1 = 1.0, 0.0, 0.0

    # If any Z value is negative, set it to 0 and renormalize the other two to sum to 1
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


def run_naive(input_prefix, out_prefix):
    run_implementation(
        compute_ibd, input_prefix, out_prefix, implementation_name="naive"
    )
