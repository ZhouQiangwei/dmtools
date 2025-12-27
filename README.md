# DMtools (dmtools)

DMtools is a high-performance toolkit and file format for DNA methylation data.
It provides a **bigWig-like, genome-indexed binary container (DM)** plus a set of CLI utilities for
fast random access, region aggregation, profiling, and **single-cell workflows** (QC, matrix, export).

> Core idea: treat methylation as a genome signal that must be queried by genomic regions efficiently,
while preserving essential methylation attributes (coverage, context, strand, optional cell ID).

---

## Highlights

- **DM format (genome-indexed, compressed, random-access)**
  - Designed for fast region queries and aggregation (similar spirit to bigWig indexing).
- **Robust conversion**
  - Convert mapping results (BAM) into DM with multi-threaded chunking (chrom/bin) and validation.
- **Bulk analysis utilities**
  - Validate files, query intervals, compute region statistics, generate profiles/summary tables.
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

### 1) Convert BAM → DM (bulk)

```bash
# Reference FASTA must match your alignment reference (same contig names).
./dmtools bam2dm \
  -g ref.fa \
  -b sample.sorted.bam \
  -o sample.dm \
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
  -b regions.bed \
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
  -o sc.sample.dm \
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
  -i sc.sample.dm \
  --idmap sc.sample.dm.idmap.tsv \
  -o sc.qc.tsv

# build sparse cell×region matrix
./dmtools sc-matrix \
  -i sc.sample.dm \
  --idmap sc.sample.dm.idmap.tsv \
  --bed regions.bed \
  --out sc_matrix \
  --sparse

# optional: empirical Bayes shrinkage (outputs eb_mean + var/CI files)
# ./dmtools sc-matrix -i sc.sample.dm --idmap sc.sample.dm.idmap.tsv --bed regions.bed --out sc_matrix --sparse --eb
```

Export to h5ad (optional Python dependency):

```bash
./dmtools sc-export \
  --prefix sc_matrix \
  --to h5ad \
  --output sc_matrix.h5ad
```

If Python/anndata is not available, sc-export should still produce the MTX bundle
(matrix.mtx, barcodes.tsv, features.tsv, plus optional QC tables) and print instructions.

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
