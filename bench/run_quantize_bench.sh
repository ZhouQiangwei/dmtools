#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bench/run_quantize_bench.sh --default-dm <path> --quantized-dm <path> --regions-bed <bed> [options]

Options:
  --output-dir <dir>   Output directory (default: bench/results/quantize-<timestamp>)
  --repeats <N>        Repeats for region query (default: 3)
  --sc-matrix          Also benchmark sc-matrix export (requires single-cell DM with idmap)
USAGE
}

DEFAULT_DM=""
QUANT_DM=""
REGIONS_BED=""
OUT_DIR=""
REPEATS=3
SC_MATRIX=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --default-dm)
      DEFAULT_DM="$2"
      shift 2
      ;;
    --quantized-dm)
      QUANT_DM="$2"
      shift 2
      ;;
    --regions-bed)
      REGIONS_BED="$2"
      shift 2
      ;;
    --output-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --repeats)
      REPEATS="$2"
      shift 2
      ;;
    --sc-matrix)
      SC_MATRIX=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
 done

if [[ -z "$DEFAULT_DM" || -z "$QUANT_DM" || -z "$REGIONS_BED" ]]; then
  usage
  exit 1
fi

RUN_ID="quantize-$(date +%Y%m%d-%H%M%S)"
if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="bench/results/$RUN_ID"
fi
mkdir -p "$OUT_DIR"

CONFIG_PATH="$OUT_DIR/config.yaml"
cat > "$CONFIG_PATH" <<CONFIG
run_id: "$RUN_ID"
paths:
  default_dm: "$DEFAULT_DM"
  quantized_dm: "$QUANT_DM"
  regions_bed: "$REGIONS_BED"
  output_dir: "$OUT_DIR/outputs"
benchmarks:
  storage:
    enabled: true
    datasets: [default_dm, quantized_dm]
  region_query:
    enabled: true
    datasets: [default_dm, quantized_dm]
    repeats: $REPEATS
    regions: ["bed"]
    command: "./dmtools regionstats -i {dm} -g {regions_bed} -o {output}"
  sc_matrix:
    enabled: $( [[ $SC_MATRIX -eq 1 ]] && echo true || echo false )
    datasets: [default_dm, quantized_dm]
    repeats: 1
    command: "./dmtools sc-matrix -i {dm} -o {output_prefix} --bed {regions_bed} --sparse"
CONFIG

python3 bench/run_benchmarks.py --config "$CONFIG_PATH"

echo "Benchmark outputs: $OUT_DIR"
