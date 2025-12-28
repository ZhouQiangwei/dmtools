# DMtools (dmtools)

DMtools is a high-performance toolkit and file format for DNA methylation data.
It provides a **bigWig-like, genome-indexed binary container (DM)** plus a set of CLI utilities for
fast random access, region aggregation, profiling, and **single-cell workflows** (QC, matrix, export).

> Core idea: treat methylation as a genome signal that must be queried by genomic regions efficiently,
while preserving essential methylation attributes (coverage, context, strand, optional cell ID).

The supported entrypoint for all core functionality is **`dmtools <command>`**
(convert, view, validate, regionstats, single-cell, DMR, export). Standalone binaries remain available
for compatibility but are not the documented interface.

---

## Highlights

- **DM format (genome-indexed, compressed, random-access)**
  - Designed for fast region queries and aggregation (similar spirit to bigWig indexing).
- **Robust conversion**
  - Convert mapping results (BAM) into DM with multi-threaded chunking (chrom/bin) and validation.
- **Bulk analysis utilities**
  - Validate files, query intervals, compute region statistics, generate profiles/summary tables.
- **DMR calling (indexed + multi-resolution)**
  - Scan DM files directly with tiled queries, merge significant tiles, then refine boundaries with finer windows and replicate-aware quasi-binomial GLM tests.
- **Single-cell support**
  - Use numeric cell/group IDs with a sidecar dictionary (`.idmap.tsv`) to scale to large cell numbers.
  - Generate sparse matrices (cell × region) and export to downstream ecosystems (e.g., h5ad).
- **Reliability**
  - Built-in validation; fail-fast on empty outputs; improved logging for merge/cleanup.

---

## Installation

### Option 1: Build from source (recommended)

```bash
git clone https://github.com/ZhouQiangwei/dmtools.git
cd dmtools

# Build bundled libraries (htslib/libBigWig) and dmtools
make libs
make

# Optional: run quick self-check
./dmtools --help
```

Note: if you build from a zip or on some filesystems, executable bits may be lost.
If you see Permission denied on htslib/version.sh, run:

```bash
chmod +x htslib/version.sh || true
make libs && make
```

---

## Quickstart

Looking for a copy/paste walkthrough with toy inputs? See [docs/quickstart.md](docs/quickstart.md); the CI
smoke test executes every command listed there.

### 1) Convert BAM → DM (bulk)

```bash
# Reference FASTA must match your alignment reference (same contig names).
./dmtools bam2dm \
  -g ref.fa \
  -b sample.sorted.bam \
  -m sample.dm \
  --threads 8 \
  --chunk-by bin \
  --validate-output
```

Validate and query:

```bash
./dmtools validate -i sample.dm

# View a small genomic region (example)
./dmtools view -i sample.dm -r chr1:100000-101000
```

Region aggregation example (BED):

```bash
./dmtools regionstats \
  -i sample.dm \
  --bed regions.bed \
  -o sample.regionstats.tsv
```

### 2) Convert BAM → DM (single-cell)

Single-cell mode requires:

- writing numeric IDs into DM records, and
- emitting a dictionary mapping id -> barcode/name in a sidecar file.

```bash
./dmtools bam2dm \
  -g ref.fa \
  -b sc.sample.sorted.bam \
  -m sc.sample.dm \
  --threads 8 \
  --chunk-by bin \
  --validate-output \
  --cell-tag CB \
  --idmap-out sc.sample.dm.idmap.tsv \
  --require-cell-tag
```

Run single-cell utilities (examples):

```bash
# per-cell QC summary
./dmtools sc-qc \
  -i sc.sample.dm -o sc.qc.tsv

# build sparse cell×region matrix
./dmtools sc-matrix \
  -i sc.sample.dm \
  -o sc_matrix \
  --bed regions.bed \
  --sparse

# optional: empirical Bayes shrinkage (outputs eb_mean + var/CI files)
# ./dmtools sc-matrix -i sc.sample.dm -o sc_matrix --bed regions.bed --sparse --eb
```

Export to h5ad (optional Python dependency):

```bash
./dmtools sc-export \
  -i sc.sample.dm \
  -o sc_matrix \
  --bed regions.bed \
  --to h5ad
```

If Python/anndata is not available, sc-export should still produce the MTX bundle
(matrix.mtx, barcodes.tsv, features.tsv, plus optional QC tables) and print instructions.

### Visualization: DM → bigWig (UCSC/IGV)

Three commands to get browser-ready tracks (example uses the tiny `test/test.dm`):

```bash
# 1) Export methylation + coverage bigWigs (ID-aware, optional coverage)
./dmtools export bigwig -i test/test.dm -o examples/trackhub_demo/out/bw --prefix demo --coverage

# 2) Build a minimal track hub (also exports bigWigs if missing)
./dmtools trackhub build --dm test/test.dm --genome demo --out-dir examples/trackhub_demo/out --with-coverage --prefix demo

# 3) Point UCSC/IGV to examples/trackhub_demo/out/hub.txt (or host that folder)
```

For merged DM files, add `--id <ID|name>` (repeatable) or `--id-list <FILE>`;
`dmtools` will split per-ID bigWigs using the provided/default idmap
(`<dm>.idmap.tsv`). A self-contained walk-through lives in
`examples/trackhub_demo/README.md`.

### Single-cell exports & reports

- Build matrices + layers: `dmtools sc-export -i sc.dm -o sc_export --bed regions.bed` (emits mC/cov/post_mean/post_var when applicable).
- Pseudobulk/cluster tracks: `dmtools sc-pseudobulk -i sc.dm -o clusters --groups groups.tsv --format bigwig --coverage`
  - Produces one bigWig (and optional coverage bigWig) per group/cluster for IGV/UCSC.
- QC report (SVG): `dmtools sc-report --qc sc_export.obs_qc.tsv --features sc_export.features.tsv --out sc_report.svg`
  - Example outputs: `docs/images/sc_qc_example.svg` (QC) and `docs/images/cluster_track_example.svg` (cluster tracks).

### DM-indexed DMR calling (coarse→refine, replicate-aware)

- Run a two-stage caller that **tiles the genome directly from DM files**, keeps FDR-stable significant tiles, then **refines boundaries with finer windows**.
- Statistics: replicate-aware **quasi-binomial logistic GLM** (weights = coverage; dispersion learned per tile/window) returning effect size, p-value, and BH-FDR q-value.
- Example:

```bash
./dmtools dmr-indexed \
  --group1 treatment.rep1.dm,treatment.rep2.dm \
  --group2 control.rep1.dm,control.rep2.dm \
  --tile-size 1000 --refine-step 250 \
  --min-delta 0.1 --coarse-q 0.1 --refine-p 0.05 \
  --out treatment_vs_control.dmr.tsv
```

Output columns: chrom/start/end, delta (mu_group1 - mu_group2), p/q, group means, dispersion (phi), methylated/coverage totals, and the number of contributing coarse tiles and refinement windows.

### Quantized mode (optional)

Use `--quantize-ml` to store values as `uint16` with scale `--quantize-scale <N>` (default `10000`).
Values are clamped to `[0,1]` and decoded as `q/scale`, so the absolute error is bounded by `1/scale`.
Smaller scales reduce file size and improve cache locality at the cost of more rounding error.

Applies to `bam2dm` and `mr2dm`, for example:

```bash
./dmtools mr2dm -C -S --Cx --quantize-ml --quantize-scale 5000 -g genome.len -m input.simpleme -o sample.quant.dm -f simpleme
```

---

## DM format overview

DM is a genome-indexed methylation container optimized for:

- random access by genomic coordinates (region queries),
- efficient region aggregation (e.g., promoters, enhancers, tiles),
- optional multiplexing by numeric ID (single-cell / group).

### ID and .idmap.tsv

For single-cell datasets, DM records store numeric uint32 IDs (e.g., 1..N),
and the mapping to cell barcodes/names is stored in:

`<file>.dm.idmap.tsv` (recommended name)

Format:

```
<id>\t<name>
1\tAAACCTGAG...
2\tAAACCTGTT...
...
```

If DM contains non-zero IDs but the idmap file is missing, single-cell commands will fail.

---

## Command overview

Run:

```bash
./dmtools --help
./dmtools <command> --help
```

Common commands (may vary slightly by branch; check --help):

- `bam2dm`: BAM → DM conversion
- `mr2dm`: text methylation ratio/count → DM (if supported)
- `validate`: DM integrity check
- `view`: query DM for a genomic region
- `regionstats`: aggregate methylation over regions (BED)
- `profile` / `bodystats` / plotting utilities (if present)
- `sc-qc`, `sc-matrix`, `sc-aggregate`, `sc-export`: single-cell workflow

---

## Troubleshooting

### Output DM is empty (or only a tiny header)

Use `--validate-output` during conversion:

dmtools will fail-fast and print why nothing was written.

Check contig naming consistency:

`ref.fa` contigs must match BAM header contigs (e.g., chr1 vs 1).

### Contig order mismatch between BAM and FASTA

dmtools uses reference contig order (FAI) for DM chrom indexing.

If your BAM contains contigs not in the reference FASTA, conversion will error (fail-fast).

### Missing .idmap.tsv in single-cell workflows

If you used `--cell-tag`/ID mode, ensure `--idmap-out` was generated and passed to sc commands.

### NFS temporary/parts directory cleanup issues

On NFS, `.nfs*` files may prevent deletion if handles are still open.

Prefer debugging with parts kept, then cleanup later.

---

## Reproducibility & Benchmarks

We recommend reporting:

- storage size under equivalent information content,
- random region query speed,
- region aggregation speed,
- single-cell matrix construction speed/memory.

See `bench/` for scripts (if included in your branch).

---

## Citation

If you use DMtools/DM in research, please cite the corresponding preprint/paper (add once available).
You may also cite the GitHub repository release you used.

---

## License

Specify your license here (e.g., MIT/BSD/GPL).
If parts depend on third-party libraries, keep their licenses in LICENSES/ as needed.

---

## Contact

Author: ZhouQiangwei
Issues and feature requests: GitHub Issues
