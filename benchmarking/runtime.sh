#!/usr/bin/env bash
# benchmark.sh — Time PLINK and python_ibd on subset and full datasets
# Keeps .log files from python_ibd runs for per-stage timing breakdown.
# Usage: bash benchmark.sh <output_directory>

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: bash benchmark.sh <output_directory>"
    exit 1
fi

BENCH_DIR="$1"
mkdir -p "$BENCH_DIR"

RESULTS="$BENCH_DIR/benchmark_results.txt"

bench() {
    local label="$1"; shift
    local start end elapsed
    
    start=$(date +%s.%N)
    
    "$@"
    
    end=$(date +%s.%N)
    
    elapsed=$(awk "BEGIN {printf \"%.3f\", $end - $start}")
    
    echo "$label: ${elapsed} s" | tee -a "$RESULTS"
}

echo "Benchmark run: $(date)" > "$RESULTS"
echo "" >> "$RESULTS"

bench "PLINK subset"       plink --bfile data/subset --genome --out "$BENCH_DIR/plink_subset"
rm -f "$BENCH_DIR/plink_subset.genome" "$BENCH_DIR/plink_subset.nosex" "$BENCH_DIR/plink_subset.log"

bench "PLINK full"         plink --bfile ~/public/ps2/ibd/ps2_ibd.lwk --genome --out "$BENCH_DIR/plink_full"
rm -f "$BENCH_DIR/plink_full.genome" "$BENCH_DIR/plink_full.nosex" "$BENCH_DIR/plink_full.log"

bench "Naive subset"       ./python_ibd --input data/subset --out "$BENCH_DIR/naive_subset" --naive
rm -f "$BENCH_DIR/naive_subset.genome"

bench "Optimized subset"   ./python_ibd --input data/subset --out "$BENCH_DIR/optimized_subset"
rm -f "$BENCH_DIR/optimized_subset.genome"

bench "Naive full"         ./python_ibd --input ~/public/ps2/ibd/ps2_ibd.lwk --out "$BENCH_DIR/naive_full" --naive
rm -f "$BENCH_DIR/naive_full.genome"

bench "Optimized full"     ./python_ibd --input ~/public/ps2/ibd/ps2_ibd.lwk --out "$BENCH_DIR/optimized_full"
rm -f "$BENCH_DIR/optimized_full.genome"

echo ""
echo "Results written to $RESULTS"
echo "Python .log files preserved in $BENCH_DIR/ for per-stage timing breakdown."