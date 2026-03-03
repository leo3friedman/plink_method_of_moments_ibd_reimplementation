# Plink Method of Moments IBD Reimplementation

## Overview

This project provides a Python implementation of the Methods of Moments algorithm for identity-by-descent (IBD) estimation found in the PLINK 1.9 toolset. The implementation relies on the statistical framework defined by [Purell et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC1950838/) and the original [PLINK 1.9 source code](https://github.com/chrchang/plink-ng). It calculates the probability of individuals sharing 0, 1, or 2 alleles via inheritance (Z0, Z1, Z2) along with the overall proportion of alleles shared (π^).

## Getting Started

1. Clone the repository using `git clone https://github.com/leo3friedman/plink_method_of_moments_ibd_reimplementation.git`.
2. Navigate into the directory with `cd plink_method_of_moments_ibd_reimplementation`.
3. Create a virtual environment via `python -m venv .env`.
4. Activate the environment using `source .env/bin/activate` for Unix/macOS or `.env\Scripts\activate` for Windows.
5. Install dependencies with `pip install -r requirements.txt`.
6. Grant execution privileges using `chmod +x python_ibd`.
7. Run a small test example with `./python_ibd --input data/subset --out python_ibd.subset`. The IBD output will be written to `python_ibd.subset.genome`.

## Usage

Run the tool using the following syntax: `./python_ibd --input data/in_prefix --out out_prefix`.

### Flags

- **`--input` (Required):** Provide the input prefix (e.g., `data/in_prefix`). The tool automatically assumes the presence of corresponding `.bed`, `.bim`, and `.fam` files in that location. See [PLINK 1.9 specifications](https://www.cog-genomics.org/plink/1.9/formats) for more information on these file types.
- **`--out` (Required):** Specify the output directory and prefix (e.g., `out_dir/out_prefix`).
- **`--naive` (Optional):** Include this flag to run the non-optimized version of the algorithm instead of the default optimized version.

### Output Files

| File Extension | Description                                        | Format                                                                                                |
| -------------- | -------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **`.genome`**  | Contains the computed IBD scores.                  | Plain text. Matches [PLINK 1.9 specification](https://www.cog-genomics.org/plink/1.9/formats#genome). |
| **`.log`**     | Contains log statements from the script execution. | Plain text.                                                                                           |

### Assumed Output Values

While the tool generates a valid `.genome` file according to the [PLINK 1.9 specification](https://www.cog-genomics.org/plink/1.9/formats#genome), it assumes specific default values for fields not calculated by this reimplementation.

| Field                 | Assumed Value |
| --------------------- | ------------- |
| `RT`                  | `'UN'`        |
| `EZ`                  | `'NA'`        |
| `PHE`                 | `-1`          |
| `DST`, `PPC`, `RATIO` | `0`           |

## Provided Datasets

The `data/` directory contains small-scale datasets and their expected PLINK `.genome` outputs for testing and validation. The full-size dataset is excluded from this repository due to size constraints, but its expected output is included as a fixture for testing.

| Dataset    | Chromosomes | Variants | Samples | Input Files                              | PLINK output          |
| :--------- | :---------- | :------- | :------ | :--------------------------------------- | :-------------------- |
| **Micro**  | 1           | 10       | 2       | `micro.bed`, `micro.bim`, `micro.fam`,   | `plink.micro.genome`  |
| **Subset** | 1           | 2000     | 10      | `subset.bed`, `subset.bim`, `subset.fam` | `plink.subset.genome` |

_Note:_ The provided datasets are subsets of the `~/public/ps2/ibd/ps2_ibd.lwk` dataset found on datahub.ucsd.edu.

## Running the Tests

The test suite verifies the accuracy of both the optimized and `--naive` implementations by comparing their output against the established PLINK 1.9 results.

- The suite actively runs both implementations against the `micro` and `subset` datasets.
- To validate the algorithm at scale without committing massive files to version control, the tests use pre-computed `.genome` fixture files to verify outputs for the full dataset.

You can execute the test suite from the root directory using `pytest` from the root directory:

`python3 -m pytest tests -v`

## Notes for Peer Reviewers

### Project Structure

```
├── python_ibd              # CLI entry point
├── requirements.txt        # Python dependencies
├── src/
│   ├── cli.py              # Argument parsing (--input, --out, --naive)
│   ├── naive.py            # Non-optimized IBD computation (loop-based)
│   ├── optimized.py        # Optimized IBD computation (vector-based)
│   └── utils.py            # Shared I/O, logging, progress, .genome writing
├── tests/
│   ├── test_naive.py       # Tests for naive implementation
│   ├── test_optimized.py   # Tests for optimized implementation
│   ├── utils.py            # Test helpers (R² calculation, .genome parsing, assertions)
│   └── fixtures/           # Pre-computed .genome files for full-dataset validation
└── data/
    ├── micro.*             # Micro dataset (2 samples, 10 variants)
    ├── subset.*            # Subset dataset (10 samples, 2000 variants)
    └── plink.*.genome      # Expected PLINK outputs for each dataset
```

### Future Work and Challenges
