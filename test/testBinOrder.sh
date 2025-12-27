#!/usr/bin/env bash
set -euo pipefail

FASTA="test/bin_order.fa"
BAM="test/bin_order.bam"
DM="test/bin_order.dm"

./test/testBinOrder "${FASTA}" "${BAM}"

./bam2dm-asan -g "${FASTA}" -b "${BAM}" --chunk-by bin --bin-size 4 -o "${DM}" --cleanup-parts

./dmtools validate -i "${DM}" --region chrA:1-10

rm -f "${FASTA}" "${FASTA}.fai" "${BAM}" "${BAM}.bai" "${DM}"
