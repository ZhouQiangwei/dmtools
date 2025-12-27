#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR=${1:-""}
if [[ -z "$RESULTS_DIR" ]]; then
  echo "Usage: bench/make_figures.sh <results_dir>" >&2
  exit 1
fi

python3 bench/make_figures.py --results-dir "$RESULTS_DIR"
