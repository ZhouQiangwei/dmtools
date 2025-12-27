# FAQ

## Chromosome naming mismatches (chr1 vs 1)

`dmtools` expects the chromosome names in your DM files, BED files, and reference to match exactly. If your reference uses `1` but your BED uses `chr1`, normalize one side (e.g., `sed 's/^chr//'`).

## TID mismatch errors

TID mismatches usually mean the contig order in your DM file does not match the order in your reference or BED. Re-generate DM files using a consistent reference and ensure all downstream inputs share the same contig ordering.

## NFS parts cleanup for chunked pipelines

Chunked BAM→DM workflows may write `.parts` intermediates. On shared or NFS filesystems, prefer explicit cleanup once you validate outputs:

```bash
find /path/to/output -name "*.parts" -delete
```

If you need to retain parts for debugging, move them to local storage before cleanup to avoid NFS latency.
