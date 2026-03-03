import itertools
import math

import numpy as np
import pandas as pd
from src.logging import StageLogger
from src.shared import (
    compute_expected_ibs,
    run_implementation,
)


def compute_ibd(genotypes: pd.DataFrame) -> np.ndarray:
    """Computes the IBD estimates (Z0, Z1, Z2) for all individuals in the genotype matrix using the method of moments approach.
    This is the naive (loop-based) implementation which is not optimized in terms of runtime.
    """

    log = StageLogger("naive.compute_ibd")

    NUM_VARIANTS = genotypes.shape[1]
    NUM_INDIVIDUALS = genotypes.shape[0]

    log.debug(
        f"Input genotype matrix has {NUM_INDIVIDUALS} individuals and {NUM_VARIANTS} variants."
    )

    # Step 1. Compute reference and alternate allele counts for each SNP.

    log.stdout(f"Stage 1/3 Computing allele frequencies...")

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

        log.update_progress(variant_id + 1, NUM_VARIANTS)

    log.finish(f"Stage 1/3 Computing allele frequencies... Done.")

    # Step 2: Compute average global expected counts of IBS states conditional on IBD states.

    log.stdout("Stage 2/3 Computing global expected IBS counts...")

    # Step 2.1: Compute expected IBS counts for polymorphic variants.
    # Note: sum_eij is equivalent to N(I=j | Z=i), (i.e. the global expected count of SNPs with IBS state j given IBD state i)

    sum_e00 = 0.0
    sum_e01 = 0.0
    sum_e02 = 0.0
    sum_e11 = 0.0
    sum_e12 = 0.0
    cnt_poly = 0

    for i, variant_stat in enumerate(variant_stats.values()):

        log.update_progress(i + 1, NUM_VARIANTS)

        X = variant_stat["ref_count"]
        Y = variant_stat["alt_count"]
        T = X + Y

        if (
            X == 0 or Y == 0 or T < 4
        ):  # PLINK 1.9 only considers polymorphic variants with > 3 observed alleles, see https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4846
            continue

        e00, e01, e02, e11, e12 = compute_expected_ibs(X, Y)
        sum_e00 += e00
        sum_e01 += e01
        sum_e02 += e02
        sum_e11 += e11
        sum_e12 += e12
        cnt_poly += 1

    # Step 2.2: Average the expected IBS counts to get the global expected counts

    avg_e00 = sum_e00 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e01 = sum_e01 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e02 = sum_e02 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e11 = sum_e11 / cnt_poly if cnt_poly > 0 else 0.0
    avg_e12 = sum_e12 / cnt_poly if cnt_poly > 0 else 0.0

    log.finish("Stage 2/3 Computing global expected IBS counts... Done.")

    # Step 3: Compute per individual pair IBD estimation.

    log.stdout("Stage 3/3 Computing pairwise IBD...")

    total_pairs = math.comb(NUM_INDIVIDUALS, 2)
    log.debug(f"Starting computing pairwise IBD for {total_pairs} pairs...")

    result = np.zeros(
        (total_pairs, 5)
    )  # columns: [individual_1_idx, individual_2_idx, Z0, Z1, Z2]

    pairs = list(itertools.combinations(range(NUM_INDIVIDUALS), 2))
    for pair_idx, (individual_1_idx, individual_2_idx) in enumerate(pairs):

        # Step 3.1: Compute IBS counts for the current pair

        ibs_0_count, ibs_1_count, ibs_2_count = 0, 0, 0

        for snp_index in range(NUM_VARIANTS):
            g1 = genotypes[individual_1_idx, snp_index]
            g2 = genotypes[individual_2_idx, snp_index]

            if math.isnan(g1) or math.isnan(
                g2
            ):  # only consider SNPs with non-missing genotypes
                continue

            diff = abs(g1 - g2)
            if diff == 2:
                ibs_0_count += 1  # opposite homozygotes: (0,2) or (2,0)
            elif diff == 0:
                ibs_2_count += 1  # same genotype: (0,0), (1,1) or (2,2)
            else:
                ibs_1_count += 1  # one allele shared: (0,1), (1,0), (1,2) or (2,1)

        S = (
            ibs_0_count + ibs_1_count + ibs_2_count
        )  # shared SNPs with non-missing genotypes

        if S == 0:  # Skip pair if no shared SNPs with non-missing genotypes
            result[pair_idx, 0] = individual_1_idx
            result[pair_idx, 1] = individual_2_idx
            continue

        # Step 3.2: Compute expected IBS count given IBD state for the current pair

        e00 = avg_e00 * S
        e01 = avg_e01 * S
        e02 = avg_e02 * S
        e11 = avg_e11 * S
        e12 = avg_e12 * S
        e22 = 1.0 * S  # All SNPs are expected to be IBS 2 given IBD state is 2

        # Step 3.3: Compute IBD estimates using Method of Moments (Purcell et al. 2007)

        z0 = ibs_0_count / e00 if e00 > 0 else 0.0
        z1 = (ibs_1_count - z0 * e01) / e11 if e11 > 0 else 0.0
        z2 = (ibs_2_count - z0 * e02 - z1 * e12) / e22 if e22 > 0 else 0.0

        # Step 3.4: Bounding IBD estimates to valid probabilities (Purcell et al. 2007)

        # Step 3.4.1: If any Z value exceeds 1, clamp it to 1 and set the other two to 0
        if z0 > 1:
            z0, z1, z2 = 1.0, 0.0, 0.0
        if z1 > 1:
            z1, z0, z2 = 1.0, 0.0, 0.0
        if z2 > 1:
            z2, z0, z1 = 1.0, 0.0, 0.0

        # Step 3.4.2: If any Z value is negative, set it to 0 and renormalize the other two to sum to 1
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

        result[pair_idx, 0] = individual_1_idx
        result[pair_idx, 1] = individual_2_idx
        result[pair_idx, 2] = z0
        result[pair_idx, 3] = z1
        result[pair_idx, 4] = z2

        log.update_progress(
            current=pair_idx + 1,
            total=total_pairs,
        )

    log.finish("Stage 3/3 Computing pairwise IBD... Done.")

    return result


def run_naive(input_prefix, out_prefix):
    run_implementation(compute_ibd, "naive", input_prefix, out_prefix)
