# DM Format Specification

DM is a bigWig-like, genome-indexed methylation container produced by dmtools. The file layout follows the bigWig/binaMeth conventions (compressed data blocks plus an R-tree index) while adding methylation-specific optional fields.

## Coordinate system

- **On-disk coordinates:** dmtools currently writes **1-based, half-open** intervals. A single cytosine at genomic position `p` is stored as `start=p`, `end=p+1`.
- **Query semantics:** interval queries are interpreted as half-open `[start, end)` ranges. When providing regions on the CLI (e.g., `chr1:1-1000000`), interpret the start as 1-based and inclusive.

## Required and optional fields

Each data record is written in the following order. Optional fields are included only if their corresponding `BM_*` flag is set in the header `version/type_mask`.

Required:

- `start` (`uint32_t`)
- `value` (`float`, methylation ratio)

Optional (in order):

- `end` (`uint32_t`, `BM_END`)
- `coverage` (`uint16_t`, `BM_COVER`)
- `strand` (`uint8_t`, `BM_STRAND`, encoded as `0` for `+`, `1` for `-`)
- `context` (`uint8_t`, `BM_CONTEXT`, encoded as `1=CG`, `2=CHG`, `3=CHH`, `0=NA`)
- `id` (`uint32_t`, `BM_ID`, numeric cell ID)

When `BM_ID` is set, dmtools expects a sidecar `<dm>.idmap.tsv` mapping numeric IDs to cell names (`id<TAB>name`).

## Data blocks and index

- Data are written in **compressed blocks** (zlib), each with a 24-byte header and `nItems` records.
- `bufSize` in the header defines the compression buffer size used for blocks (default: **32768** bytes).
- An **R-tree index** references each compressed block to support random access. The index `blockSize` is the maximum number of child entries per node (default: **256**).
- Zoom levels are supported by the format but typically unused; the zoom headers are reserved in the file header.

## Header versioning and compatibility

- The header `version/type_mask` field stores `BM_MAGIC` OR-ed with `BM_*` flags that describe record layout.
- Readers must check that `BM_MAGIC` is present before attempting to decode records.
- New optional fields must be added as new `BM_*` flags while preserving the ordering rules above.

## Write-parameter extension

If `extensionOffset` is non-zero, dmtools writes a small extension describing write parameters (to aid reproduction and debugging):

- `magic = 0x44574d50` ("DWMP")
- `version = 1`
- `size = sizeof(struct)`
- `bufSize` and `blockSize` used by the writer

Readers that do not recognize the extension can safely ignore it.

## CLI naming guidance

- The primary interface is `dmtools <mode> ...` (for example, `dmtools bam2dm`).
- Standalone binaries (e.g., `bam2dm`) mirror the same behavior and `--help` output.
- Prefer `dmtools --help` and `dmtools <mode> --help` for consistent option documentation.
