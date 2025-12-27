#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-bench/config.yaml}
python3 bench/run_benchmarks.py --config "$CONFIG"
