# Quickstart

These examples use the `dmtools` entrypoint for every operation. Run them from the repository root **after** building `./dmtools`; each block stays under eight commands.

## Bulk methylation (toy data)

```bash
DMTOOLS=$(pwd)/dmtools; workdir=$(mktemp -d); cd "$workdir"

cat <<'EOF' > bulk.treatment.simpleme
chr1    5   0.8 10  + CG
chr1    20  0.3 12  - CHG
chr1    55  0.6 15  + CG
EOF

cat <<'EOF' > bulk.control.simpleme
chr1    5   0.4 12  + CG
chr1    20  0.5 10  - CHG
chr1    55  0.2 18  + CG
EOF

printf "chr1\t1000\n" > genome.len && printf "chr1\t1\t80\ttoy_region\n" > bulk.regions.bed

$DMTOOLS mr2dm -C -S --Cx -g genome.len -m bulk.treatment.simpleme -o bulk.treatment.dm -f simpleme
$DMTOOLS mr2dm -C -S --Cx -g genome.len -m bulk.control.simpleme -o bulk.control.dm -f simpleme
$DMTOOLS validate -i bulk.treatment.dm
$DMTOOLS view -i bulk.treatment.dm -r chr1:1-80 -o bulk.view.tsv
$DMTOOLS regionstats -i bulk.treatment.dm --bed bulk.regions.bed -o bulk.regionstats.tsv
$DMTOOLS dmr-indexed --group1 bulk.treatment.dm --group2 bulk.control.dm --tile-size 50 --refine-step 25 --out bulk.dmr.tsv
```

## Single-cell methylation (toy data)

```bash
DMTOOLS=$(pwd)/dmtools; workdir=$(mktemp -d); cd "$workdir"

cat <<'EOF' > sc.simpleme
chr1    10  0.9 12  + CG
chr1    40  0.2 8   - CHG
chr1    70  0.5 6   + CHH
EOF

printf "chr1\t1000\n" > sc.genome.len && printf "chr1\t1\t50\tregionA\nchr1\t51\t100\tregionB\n" > sc.regions.bed && printf "cell1\tgroupA\n" > cell_groups.tsv

$DMTOOLS mr2dm -C -S --Cx --Id -g sc.genome.len -m sc.simpleme -o sc.dm -f simpleme --cell-name cell1 --idmap-out sc.dm.idmap.tsv
$DMTOOLS sc-qc -i sc.dm -o sc.qc.tsv
$DMTOOLS sc-matrix -i sc.dm -o sc_matrix --bed sc.regions.bed --sparse
$DMTOOLS sc-aggregate -i sc.dm -o sc_agg --groups cell_groups.tsv --bed sc.regions.bed --dense
$DMTOOLS sc-export -i sc.dm -o sc_export --bed sc.regions.bed  # emits mC/cov/post_mean/post_var when available
$DMTOOLS sc-pseudobulk -i sc.dm -o clusters --groups cell_groups.tsv --format bigwig --coverage
$DMTOOLS sc-report --qc sc_export.obs_qc.tsv --features sc_export.features.tsv --out sc_report.svg
```

## Quantized mode (optional)

For smaller files, add `--quantize-ml --quantize-scale <N>` when writing (`bam2dm` or `mr2dm`).
Values are stored as `uint16` and decoded as `value = q/scale`, giving an absolute error ≤ `1/scale`.
Example:

```bash
$DMTOOLS mr2dm -C -S --Cx --quantize-ml --quantize-scale 5000 -g genome.len -m toy.simpleme -o toy.quant.dm -f simpleme
```
