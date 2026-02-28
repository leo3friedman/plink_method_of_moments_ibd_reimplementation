import argparse


def parse_args():
    """Parses command-line arguments for the IBD calculation tool."""
    parser = argparse.ArgumentParser(
        description="Compute pairwise Identity-by-Descent (IBD) from PLINK files (.bed, .fam, .bim).",
        epilog="Example usage: ./python_ibd --input data/sample_data --out tmp/sample.ibd",
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Prefix for the input files (assumes .bed, .bim, .fam exist). Example: data/in_prefix",
    )

    parser.add_argument(
        "--out",
        required=True,
        type=str,
        help="Prefix for the output files (.genome and .log will be generated). Example: out/out_prefix",
    )

    parser.add_argument(
        "--naive",
        action="store_true",
        help="Run the naive implementation. If omitted, defaults to the optimized implementation.",
    )

    return parser.parse_args()
