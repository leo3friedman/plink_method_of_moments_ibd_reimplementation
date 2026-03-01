from src.naive import run_naive
from tests.utils import (
    assert_implementation_matches_plink,
    assert_log_file_exists,
    get_tmp_output_prefix,
)


def test_naive_matches_plink_on_subset():
    input_prefix = "data/subset"
    ground_truth = "data/plink.subset.genome"
    out_prefix = get_tmp_output_prefix("naive_test_subset")
    assert_implementation_matches_plink(
        input_prefix, out_prefix, ground_truth, run_naive
    )
    assert_log_file_exists(out_prefix)


def test_naive_matches_plink_on_micro():
    input_prefix = "data/micro"
    ground_truth = "data/plink.micro.genome"
    out_prefix = get_tmp_output_prefix("naive_test_micro")
    assert_implementation_matches_plink(
        input_prefix, out_prefix, ground_truth, run_naive
    )
    assert_log_file_exists(out_prefix)
