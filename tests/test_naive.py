from src.naive import run_naive
from tests.shared import (
    assert_log_file_exists,
    get_tmp_output_prefix,
    assert_genome_data,
)


def test_naive_matches_plink_on_micro():
    input_prefix = "data/micro"
    out_prefix = get_tmp_output_prefix("naive_test_micro")
    run_naive(input_prefix, out_prefix)
    ground_truth_filepath = "data/plink.micro.genome"
    target_genome_filepath = f"{out_prefix}.genome"

    assert_genome_data(
        truth_filepath=ground_truth_filepath,
        target_filepath=target_genome_filepath,
        r_squared_threshold=0.99,
    )
    assert_log_file_exists(out_prefix)


def test_naive_matches_plink_on_subset():
    input_prefix = "data/subset"
    out_prefix = get_tmp_output_prefix("naive_test_subset")
    run_naive(input_prefix, out_prefix)
    ground_truth_filepath = "data/plink.subset.genome"
    target_genome_filepath = f"{out_prefix}.genome"
    assert_genome_data(
        truth_filepath=ground_truth_filepath,
        target_filepath=target_genome_filepath,
        r_squared_threshold=0.99,
    )
    assert_log_file_exists(out_prefix)


def test_naive_matches_plink_on_full():
    ground_truth = "tests/fixtures/plink.full.genome"
    target = "tests/fixtures/python_ibd.naive.full.genome"  # use a fixed output prefix for the full test since it takes a long time to run

    assert_genome_data(
        truth_filepath=ground_truth,
        target_filepath=target,
        r_squared_threshold=0.999,  # use a higher threshold for the full dataset since it has more pairs and thus less noise
    )
