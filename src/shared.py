import os
import pandas as pd
from bed_reader import open_bed
from src.logging import StageLogger


def compute_expected_ibs(X, Y):
    """Compute per-variant expected IBS contributions from ref count X and alt count Y.

    Callers must filter out variants where X == 0, Y == 0, or X + Y < 4
    before calling this function to avoid division by zero.

    Works with both scalar values (loop-based) and numpy arrays (vectorized).

    See Purcell et al. 2007 Table 1 and PLINK 1.9:
    https://github.com/chrchang/plink-ng/blob/c785858ab8ebfd62fe8367d9a878323607086fde/1.9/plink_calc.c#L4849-L4866
    """

    T = X + Y
    p = X / T
    q = Y / T
    dpp_sq = p**2
    dqq_sq = q**2
    dxx1 = (X - 1) / X
    dxx2 = (X - 1) * (X - 2) / (X**2)
    dyy1 = (Y - 1) / Y
    dyy2 = (Y - 1) * (Y - 2) / (Y**2)
    num_allelesf2 = T * T / ((T - 1) * (T - 2))
    num_allelesf3 = num_allelesf2 * T / (T - 3)

    e00 = 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3
    e01 = 4 * p * q * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2)
    e02 = num_allelesf3 * (
        dqq_sq * dqq_sq * dyy2 * (Y - 3) / Y
        + dpp_sq * dpp_sq * dxx2 * (X - 3) / X
        + 4 * dpp_sq * dqq_sq * dxx1 * dyy1
    )
    e11 = 2 * p * q * num_allelesf2 * (p * dxx1 + q * dyy1)
    e12 = num_allelesf2 * (
        dpp_sq * p * dxx2 + dqq_sq * q * dyy2 + dpp_sq * q * dxx1 + p * dqq_sq * dyy1
    )
    return e00, e01, e02, e11, e12


def run_implementation(
    ibd_fn: callable,
    implementation_name: str,
    input_prefix: str,
    out_prefix: str,
):

    StageLogger.setup(out_prefix)

    log = StageLogger("utils.run_implementation")
    log.debug(f"input prefix: '{input_prefix}', output prefix: '{out_prefix}'")

    log.stdout("Reading input file...")

    bed = open_bed(f"{input_prefix}.bed")
    individual_ids = bed.iid
    family_ids = bed.fid
    genotypes = bed.read()

    log.finish("Reading input file... Done.")

    log.debug(
        f"Input files contain {genotypes.shape[0]} individuals, {genotypes.shape[1]} variants."
    )

    log.stdout(f"Computing pairwise IBD using {implementation_name} implementation...")

    ibd_results = ibd_fn(genotypes)

    idx1 = ibd_results[:, 0].astype(int)
    idx2 = ibd_results[:, 1].astype(int)
    z0 = ibd_results[:, 2]
    z1 = ibd_results[:, 3]
    z2 = ibd_results[:, 4]

    result = pd.DataFrame(
        {
            "FID1": family_ids[idx1],
            "IID1": individual_ids[idx1],
            "FID2": family_ids[idx2],
            "IID2": individual_ids[idx2],
            "RT": "UN",
            "EZ": "NA",
            "Z0": z0,
            "Z1": z1,
            "Z2": z2,
            "PI_HAT": z1 / 2 + z2,
            "PHE": -1,
            "DST": 0,
            "PPC": 0,
            "RATIO": 0,
        }
    )

    log.stdout("Writing output files...")

    # make file if does not exist, otherwise overwrite
    filepath = f"{out_prefix}.genome"

    dir_name = os.path.dirname(filepath)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    with open(filepath, "w") as f:
        f.write(result.to_string(justify="right", index=False))

    log.finish(f"Output written to {out_prefix}.genome and {out_prefix}.log")
