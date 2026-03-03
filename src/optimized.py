import pandas as pd
import numpy as np
import math

from src.logging import StageLogger

from src.shared import (
    compute_expected_ibs,
    run_implementation,
)


def compute_ibd(genotypes: pd.DataFrame) -> pd.DataFrame:
    """Computes the IBD estimates (Z0, Z1, Z2) for all individuals in the genotype matrix using the method of moments approach.
    This is the optimized (vector-based) implementation which should outperform naive.compute_ibd on large datasets in terms of runtime.
    """

    log = StageLogger("optimized.compute_ibd")

    NUM_VARIANTS = genotypes.shape[1]
    NUM_INDIVIDUALS = genotypes.shape[0]

    log.debug(
        f"Input genotype matrix has {NUM_INDIVIDUALS} individuals and {NUM_VARIANTS} variants."
    )

    # Step 1: Compute reference and alternate allele counts for each SNP.

    log.stdout(f"Stage 1/3 Computing allele frequencies...")

    non_missing = ~np.isnan(genotypes)
    alt_counts = np.nansum(genotypes, axis=0)
    ref_counts = 2 * non_missing.sum(axis=0) - alt_counts

    log.finish(f"Stage 1/3 Computing allele frequencies... Done.")

    # Step 2: Compute average global expected counts of IBS states conditional on IBD states.

    log.stdout("Stage 2/3 Computing global expected IBS counts...")

    # Step 2.1: Compute expected IBS counts for polymorphic variants.

    T = ref_counts + alt_counts
    poly = (
        (ref_counts > 0) & (alt_counts > 0) & (T > 3)
    )  # PLINK 1.9 only considers polymorphic variants with > 3 observed alleles, see https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4846

    e00_v, e01_v, e02_v, e11_v, e12_v = compute_expected_ibs(
        ref_counts[poly], alt_counts[poly]
    )

    # Step 2.2: Average the expected IBS counts to get the global expected average per SNP
    # Note: avg_eij is the global expectation of seeing IBS state j given IBD state i for any arbitrary SNP

    cnt_poly = poly.sum()
    avg_e00 = e00_v.sum() / cnt_poly if cnt_poly > 0 else 0.0
    avg_e01 = e01_v.sum() / cnt_poly if cnt_poly > 0 else 0.0
    avg_e02 = e02_v.sum() / cnt_poly if cnt_poly > 0 else 0.0
    avg_e11 = e11_v.sum() / cnt_poly if cnt_poly > 0 else 0.0
    avg_e12 = e12_v.sum() / cnt_poly if cnt_poly > 0 else 0.0

    log.finish("Stage 2/3 Computing global expected IBS counts... Done.")

    # Step 3: Compute per-individual pair IBD estimation.

    log.stdout("Stage 3/3 Computing pairwise IBD...")

    total_pairs = math.comb(NUM_INDIVIDUALS, 2)
    log.debug(f"Computing pairwise IBD for {total_pairs} pairs.")

    # Step 3.1: Create IBS state indicators for each IBS state
    # Note: isk[i,m] = 1 if individual i has genotype k at SNP m

    is0 = (genotypes == 0).astype(np.float32)  # Shape: (n_individuals, n.variants)
    is1 = (genotypes == 1).astype(np.float32)
    is2 = (genotypes == 2).astype(np.float32)

    # Step 3.2: Compute number of non-missing SNPs per individual pair
    # Note: S[i,j] = count of SNPs with non-missing genotypes for individuals i and j.

    valid = (~np.isnan(genotypes)).astype(np.float32)
    S = valid @ valid.T  # Shape: (n_individuals, n_individuals)

    # Step 3.3: Compute IBS counts for all each genome
    # Note: ibs_k_counts[i, j] = count of IBS k SNPs between individuals i and j.
    # Note: (isk @ isk.T)[i, j] is the count of SNPs where individuals i and j both have genotype k

    ibs_0_counts = is0 @ is2.T + is2 @ is0.T  # Shape: (n_individuals, n_individuals)
    ibs_2_counts = is0 @ is0.T + is1 @ is1.T + is2 @ is2.T
    ibs_1_counts = S - ibs_0_counts - ibs_2_counts  # all remaining SNPs

    # Step 3.4: Restrict to upper triangle pairs, ibs_k[i] = counts of SNPs with IBS k for pair i (i.e. individuals (ind_i[i], ind_j[i]))
    # Note: this is equivalent to itertools.combinations. The ibs_k_counts (n_individuals, n_individuals) matrices represents all permutations of pairs, but we only care about combinations (i.e. (0,1) is the same as (1,0)).

    idx_i, idx_j = np.triu_indices(NUM_INDIVIDUALS, k=1)
    ibs_0 = ibs_0_counts[idx_i, idx_j]  # Shape: (n_pairs, )
    ibs_1 = ibs_1_counts[idx_i, idx_j]
    ibs_2 = ibs_2_counts[idx_i, idx_j]
    S = S[idx_i, idx_j]

    # Step 3.5: Compute per-pair expected IBS counts given IBD state
    # Note: eij[k] = expected count of IBS j SNPs for pair k given IBD state i

    e00 = avg_e00 * S  # Shape: (n_pairs, )
    e01 = avg_e01 * S
    e02 = avg_e02 * S
    e11 = avg_e11 * S
    e12 = avg_e12 * S
    e22 = 1.0 * S  # All SNPs are expected to be IBS 2 given IBD state is 2

    # Step 3.6: Compute IBD estimates using Method of Moments (Purcell et al. 2007)

    z0 = np.where(e00 > 0, ibs_0 / e00, 0.0)  # Shape: (n_pairs, )
    z1 = np.where(e11 > 0, (ibs_1 - z0 * e01) / e11, 0.0)
    z2 = np.where(e22 > 0, (ibs_2 - z0 * e02 - z1 * e12) / e22, 0.0)

    log.debug("Finished computing raw IBD estimates, applying bounding procedure.")

    # Step 3.7: Bounding IBD estimates to valid probabilities (Purcell et al. 2007)

    # Step 3.7.1: If any Z value exceeds 1, clamp it to 1 and set the other two to 0

    z1 = np.where(z0 > 1, 0.0, z1)
    z2 = np.where(z0 > 1, 0.0, z2)
    z0 = np.where(z0 > 1, 1.0, z0)

    z0 = np.where(z1 > 1, 0.0, z0)
    z2 = np.where(z1 > 1, 0.0, z2)
    z1 = np.where(z1 > 1, 1.0, z1)

    z0 = np.where(z2 > 1, 0.0, z0)
    z1 = np.where(z2 > 1, 0.0, z1)
    z2 = np.where(z2 > 1, 1.0, z2)

    # Step 3.7.2: If any Z value is negative, set it to 0 and renormalize the other two to sum to 1

    mask = z0 < 0
    z0 = np.where(mask, 0.0, z0)
    s = z1 + z2
    z1 = np.where(mask & (s > 0), z1 / np.where(s > 0, s, 1.0), z1)
    z2 = np.where(mask & (s > 0), z2 / np.where(s > 0, s, 1.0), z2)

    mask = z1 < 0
    z1 = np.where(mask, 0.0, z1)
    s = z0 + z2
    z0 = np.where(mask & (s > 0), z0 / np.where(s > 0, s, 1.0), z0)
    z2 = np.where(mask & (s > 0), z2 / np.where(s > 0, s, 1.0), z2)

    mask = z2 < 0
    z2 = np.where(mask, 0.0, z2)
    s = z0 + z1
    z0 = np.where(mask & (s > 0), z0 / np.where(s > 0, s, 1.0), z0)
    z1 = np.where(mask & (s > 0), z1 / np.where(s > 0, s, 1.0), z1)

    log.finish("Stage 3/3 Computing pairwise IBD... Done.")

    result = np.zeros((total_pairs, 5))
    result[:, 0] = idx_i  # individual 1 index
    result[:, 1] = idx_j  # individual 2 index
    result[:, 2] = z0
    result[:, 3] = z1
    result[:, 4] = z2

    return result


def run_optimized(input_prefix, out_prefix):
    run_implementation(compute_ibd, "optimized", input_prefix, out_prefix)
