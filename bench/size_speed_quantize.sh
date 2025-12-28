#!/usr/bin/env bash
set -euo pipefail

# Minimal quantization size/speed comparison using toy data.
# Generates float32 and quantized DM files, prints size and timing for region query + regionstats.

DMTOOLS=${DMTOOLS:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/dmtools"}
SCALE=${SCALE:-10000}
REPEATS=${REPEATS:-3}

if [[ ! -x "${DMTOOLS}" ]]; then
  echo "dmtools binary not found at ${DMTOOLS}" >&2
  exit 1
fi

workdir=$(mktemp -d)
if [[ -z "${KEEP_WORKDIR:-}" ]]; then
  trap 'rm -rf "$workdir"' EXIT
else
  echo "[bench] KEEP_WORKDIR is set; artifacts remain in $workdir"
fi
cd "$workdir"

cat > toy.simpleme <<'EOF'
chr1    10  0.10 10  + CG
chr1    50  0.35 20  - CHH
chr1    80  0.80 30  + CG
EOF
printf "chr1\t200\n" > toy.genome.len
printf "chr1\t1\t100\n" > toy.regions.bed

echo "[bench] writing float32 DM"
"${DMTOOLS}" mr2dm -C -S --Cx -g toy.genome.len -m toy.simpleme -o float.dm -f simpleme
echo "[bench] writing quantized DM (scale=${SCALE})"
"${DMTOOLS}" mr2dm -C -S --Cx -g toy.genome.len -m toy.simpleme -o quant.dm -f simpleme --quantize-ml --quantize-scale "${SCALE}"

SIZE_FLOAT=$(python3 - <<'PY'
import os
print(os.path.getsize("float.dm"))
PY
)
SIZE_QUANT=$(python3 - <<'PY'
import os
print(os.path.getsize("quant.dm"))
PY
)

measure() {
  python3 - "$@" <<'PY'
import os, subprocess, sys, time
repeats = int(os.environ.get("REPEATS", "3"))
cmd = sys.argv[1:]
times = []
for _ in range(repeats):
    t0 = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    times.append(time.perf_counter() - t0)
print(f"{sum(times)/len(times):.4f}")
PY
}

echo "[bench] timing region query (${REPEATS} repeats)"
RQ_FLOAT=$(measure "${DMTOOLS}" view -i float.dm -r chr1:1-120 -o /dev/null)
RQ_QUANT=$(measure "${DMTOOLS}" view -i quant.dm -r chr1:1-120 -o /dev/null)

echo "[bench] timing regionstats (${REPEATS} repeats)"
RS_FLOAT=$(measure "${DMTOOLS}" regionstats -i float.dm --bed toy.regions.bed -o /dev/null)
RS_QUANT=$(measure "${DMTOOLS}" regionstats -i quant.dm --bed toy.regions.bed -o /dev/null)

cat <<REPORT
Quantization comparison (toy data)
----------------------------------
Value scale: ${SCALE}
File size (bytes): float=${SIZE_FLOAT} quant=${SIZE_QUANT}
Region query time (s): float=${RQ_FLOAT} quant=${RQ_QUANT}
Regionstats time (s): float=${RS_FLOAT} quant=${RS_QUANT}
DM files: ${workdir}/float.dm , ${workdir}/quant.dm
Cleanup: $( [[ -z "${KEEP_WORKDIR:-}" ]] && echo "temporary dir auto-removed" || echo "kept (KEEP_WORKDIR set)" )
REPORT
