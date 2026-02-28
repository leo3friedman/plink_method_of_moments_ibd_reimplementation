import tempfile
from src.naive import run_naive
from tests.utils import validate_genome_data
import os


def test_naive_matches_plink():
    input_prefix = "data/subset"
    ground_truth = "data/plink.subset.genome"

    # create a temporary output file for the naive implementation
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_output_prefix = f"{tmpdir}/naive_test"

    run_naive(input_prefix, tmp_output_prefix)

    validate_genome_data(ground_truth, f"{tmp_output_prefix}.genome")

    # validate log file exists
    assert os.path.exists(
        f"{tmp_output_prefix}.log"
    ), "Expected log file to be generated but it does not exist"
