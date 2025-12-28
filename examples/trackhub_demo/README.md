# Track hub demo (bigWig export)

This folder shows a minimal, reproducible path from a DM file to UCSC/IGV-ready
bigWig tracks and a track hub.

Prerequisites:

- `dmtools` built in the repository root (`./dmtools`).

Quick demo (uses the tiny `test/test.dm`):

```bash
DMTOOLS=./dmtools
DEMO_OUT=examples/trackhub_demo/out

# 1) Export methylation + coverage bigWigs (optional when using step 2)
$DMTOOLS export bigwig -i test/test.dm -o "$DEMO_OUT/bw" --prefix demo --coverage

# 2) Build a minimal hub (also exports bigWigs to $DEMO_OUT/bw if missing)
$DMTOOLS trackhub build --dm test/test.dm --genome demo --out-dir "$DEMO_OUT" --with-coverage --prefix demo

# 3) Inspect hub.txt / genomes.txt / trackDb.txt and bigWigs under $DEMO_OUT
ls "$DEMO_OUT"
```

Upload the `out/` directory (or serve it over HTTP) and point UCSC/IGV to
`hub.txt`. Replace `test/test.dm` with your own DM (or merged DM + `--id ...`)
to generate per-sample/per-ID tracks.
