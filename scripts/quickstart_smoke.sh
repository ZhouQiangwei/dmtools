#!/usr/bin/env bash
set -euo pipefail

# Runs the exact toy commands documented in docs/quickstart.md.

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DMTOOLS="${REPO_ROOT}/dmtools"

run_block() {
  local label="$1"
  shift
  local dir
  dir=$(mktemp -d)
  echo "[quickstart] ${label} workspace: ${dir}"
  pushd "$dir" >/dev/null
  "$@"
  popd >/dev/null
  rm -rf "$dir"
}

run_block bulk bash -c '
  cat <<'"'"'EOF'"'"' > bulk.simpleme
chr1    5   0.8 10  + CG
chr1    20  0.3 12  - CHG
chr1    55  0.6 15  + CG
EOF

  printf "chr1\t1000\n" > genome.len && printf "chr1\t1\t50\ttoy_region\n" > bulk.regions.bed && printf "chr1\t1\t50\t2\t5\t3\t6\n" > bulk.dmr_counts.tsv

  "'"${DMTOOLS}"'" mr2dm -C -S --Cx -g genome.len -m bulk.simpleme -o bulk.dm -f simpleme
  "'"${DMTOOLS}"'" validate -i bulk.dm
  "'"${DMTOOLS}"'" view -i bulk.dm -r chr1:1-60 -o bulk.view.tsv
  "'"${DMTOOLS}"'" regionstats -i bulk.dm --bed bulk.regions.bed -o bulk.regionstats.tsv
  "'"${DMTOOLS}"'" dmr-bb --input bulk.dmr_counts.tsv --out bulk.dmr.tsv
'

run_block single-cell bash -c '
  cat <<'"'"'EOF'"'"' > sc.simpleme
chr1    10  0.9 12  + CG
chr1    40  0.2 8   - CHG
chr1    70  0.5 6   + CHH
EOF

  printf "chr1\t1000\n" > sc.genome.len && printf "chr1\t1\t50\tregionA\nchr1\t51\t100\tregionB\n" > sc.regions.bed && printf "cell1\tgroupA\n" > cell_groups.tsv

  "'"${DMTOOLS}"'" mr2dm -C -S --Cx --Id -g sc.genome.len -m sc.simpleme -o sc.dm -f simpleme --cell-name cell1 --idmap-out sc.dm.idmap.tsv
  "'"${DMTOOLS}"'" sc-qc -i sc.dm -o sc.qc.tsv
  "'"${DMTOOLS}"'" sc-matrix -i sc.dm -o sc_matrix --bed sc.regions.bed --sparse
  "'"${DMTOOLS}"'" sc-aggregate -i sc.dm -o sc_agg --groups cell_groups.tsv --bed sc.regions.bed --dense
  "'"${DMTOOLS}"'" sc-export -i sc.dm -o sc_export --bed sc.regions.bed
'
