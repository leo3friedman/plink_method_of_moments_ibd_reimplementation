from bed_reader import open_bed
import itertools
from src.utils import GENOME_COLUMNS, write_genome_file
import pandas as pd


def run_naive(input_prefix, out_prefix):
    bed = open_bed(f"{input_prefix}.bed")
    individual_ids = bed.iid
    family_ids = bed.fid

    sample_ids = list(zip(individual_ids, family_ids))
    pairwise_combinations = list(itertools.combinations(sample_ids, 2))
    result = pd.DataFrame(
        columns=GENOME_COLUMNS,
        data=[
            {
                "FID1": pair[0][1],
                "IID1": pair[0][0],
                "FID2": pair[1][1],
                "IID2": pair[1][0],
                "RT": "UN",
                "EZ": "NA",
                "Z0": 0,
                "Z1": 0,
                "Z2": 0,
                "PI_HAT": 0,
                "PHE": -1,
                "DST": 0,
                "PPC": 0,
                "RATIO": 0,
            }
            for pair in pairwise_combinations
        ],
    )

    output_filepath = f"{out_prefix}.genome"
    write_genome_file(result, output_filepath)
