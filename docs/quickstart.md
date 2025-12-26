# Quickstart

```bash
# Build
make -j$(nproc)

# Basic view
./dmtools view -i input.dm -r chr1:100000-101000 -o /dev/stdout

# Region aggregation
./dmtools regionstats -i input.dm -g regions.bed -o region_stats.tsv

# Single-cell matrix (sparse bundle)
./dmtools sc-matrix -i input.dm -o sc_matrix --bed regions.bed --sparse
```
