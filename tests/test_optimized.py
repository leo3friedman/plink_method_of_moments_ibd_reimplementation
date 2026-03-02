from src.optimized import run_optimized
from tests.utils import (
    assert_log_file_exists,
    get_tmp_output_prefix,
    assert_genome_data,
)


def test_optimized_matches_plink_on_micro():
    input_prefix = "data/micro"
    out_prefix = get_tmp_output_prefix("optimized_test_micro")
    run_optimized(input_prefix, out_prefix)
    ground_truth_filepath = "data/plink.micro.genome"
    target_genome_filepath = f"{out_prefix}.genome"

    assert_genome_data(
        truth_filepath=ground_truth_filepath,
        target_filepath=target_genome_filepath,
        r_squared_threshold=0.99,
    )
    assert_log_file_exists(out_prefix)
