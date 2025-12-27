<p align="left">
  <a href="https://dmtools-docs.rtfd.io" target="_blank"><img src="https://img.shields.io/badge/docs-8A2BE2?link=https%3A%2F%2Fdmtools-docs.rtfd.io"></a>
  <a href="https://github.com/ZhouQiangwei/dmtools/releases"><img src="https://img.shields.io/github/v/release/ZhouQiangwei/dmtools"></a>
  <img src="https://img.shields.io/github/license/ZhouQiangwei/dmtools">
  <img src="https://img.shields.io/github/stars/ZhouQiangwei/dmtools">
  <a href="https://github.com/ZhouQiangwei/dmtools/issues"><img src="https://img.shields.io/github/issues/ZhouQiangwei/dmtools?color=red"></a>
  <img src="https://img.shields.io/github/forks/ZhouQiangwei/dmtools?color=teal">

  <img src="https://img.shields.io/github/downloads/ZhouQiangwei/dmtools/total">
  <img src="https://img.shields.io/github/watchers/ZhouQiangwei/dmtools">
  <a href='https://github.com/ZhouQiangwei/dmtools'><img alt='Clones' src='https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://gist.githubusercontent.com/ZhouQiangwei/6f7019fe4f7d99f0eb0af4dd610a51db/raw/clone.json&logo=github'></a>
  <a href="https://dmtools-docs.rtfd.io"><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FZhouQiangwei%2Fdmtools&count_bg=%2379C83D&title_bg=%23555555&icon=coffeescript.svg&icon_color=%23E7E7E7&title=views&edge_flat=false"/></a>
</p>


# DM format and DMtools

DM is a bigWig-like, genome-indexed methylation container designed for compact storage and fast random access. DMtools is a C/C++ toolkit for generating, validating, and analyzing DM files.

DM format stores methylation ratios with optional coverage, context, strand, and ID fields (for single-cell datasets).

DMtools provides file format conversion, quick view, methylation level calculation, differential methylation analysis, and other DM-specific utilities.

## Documentation

- [Install](docs/install.md)
- [Quickstart](docs/quickstart.md)
- [FAQ](docs/faq.md)
- [Changelog](CHANGELOG.md)

For more information, please see https://dmtools-docs.rtfd.io

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://dmtools-docs.readthedocs.io/en/latest/function/bam2dm.html


## Installation and build dependencies

The native binaries are built with GNU make and require the following system packages:

- A C/C++ toolchain (GCC 5+ recommended; GCC 4.x works via the `-std=gnu++0x` fallback in the Makefile), `make`, and standard build tools
- Compression and network dependencies: `zlib`/`libz-dev`, `libbz2-dev`, `liblzma-dev`, `libcurl4-openssl-dev`
- Optional: GSL (`libgsl-dev`). When present, `dmDMR` is built; otherwise it is skipped automatically and the build prints a short notice.

### Build without root (recommended)

All targets can be installed under a user prefix without sudo:

```
export PREFIX="$HOME/.local"
mkdir -p "$PREFIX"
make clean
make libs
make
make install PREFIX="$PREFIX"
export PATH="$PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$LD_LIBRARY_PATH"
```

Verify the toolchain:

```
dmtools bam2dm --help
```

### Conda-based toolchain (no sudo)

If your cluster lacks the required system libraries, a Conda environment can provide them without root:

```
conda create -n dmtools-env -c conda-forge gcc gxx zlib bzip2 xz libcurl make
conda activate dmtools-env
make clean
make libs
make
```

### Optional system package installation

If you do have root on Debian/Ubuntu, the following installs the prerequisites. This is optional because the above user-local builds avoid sudo:

```
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
```

To include `dmDMR`, add `libgsl-dev` and rebuild with `make WITH_GSL=1`.

The build outputs the `bam2dm` executable from `bam2dm.cpp`; `dmtools bam2dm` dispatches to that binary. Use `dmtools --help` and `dmtools <mode> --help` for consistent CLI documentation; standalone binaries mirror the same options.
When you run `dmtools bam2dm`, the driver shells out to the compiled `bam2dm` in the same directory. Use `-p N` **or** `--threads N` to control how many chromosome-parallel workers the wrapper launches (default is 1). If you pass `-p/--threads N` with `N>1`, `dmtools` fans out by chromosome (`--chrom <name>`) and then merges the per-chromosome `-m <outfile>.<chrom>` shards back together with `dmtools merge`, matching the legacy behavior of the chromosome-parallel wrapper.


| [ dmtools ]         | process with dm file                                                    |
|---------------------|--------------------------------------------------------------------------|
|                     | dmtools <mode> [opnions]                                                 |
| [mode]              | index align mr2dm view overlap regionstats bodystats profile chromstats |
| index               | build index for genome                                                  |
| align               | alignment wrapper (if built; see `dmtools --help`)                       |
| mr2dm               | convert txt meth file to dm format                                      |
| view                | dm format to txt meth                                                   |
| ebsrate             | estimate bisulfite conversion rate                                      |
| viewheader          | view header of dm file                                                  |
| overlap             | overlap cytosine site with more than two dm files                       |
| regionstats         | calculate DNA methylation level of per region                            |
| bodystats           | calculate DNA methylation level of body, upstream and downstream.        |
| profile             | calculate DNA methylation profile                                        |
| chromstats          | calculate DNA methylation level across chromosome                        |
| chrmeth             | calculate DNA methylation level of all chromosomes                       |
| addzm               | add or change zoom levels for dm format, need for browser visulization   |
| stats               | coverage and methylation level distribution of data                      |
| -h|--help                                                                                      |

| dmDMR               | detect DMC/DMR based on dm file                                         |
|---------------------|--------------------------------------------------------------------------|
| -p                  | output file prefix                                                       |
| -1                  | sample1 methy dm files, sperate by comma.                               |
| -2                  | sample2 methy dm files, sperate by comma.                               |
| --mindmc            | min dmc sites in dmr region. [default : 4]                               |
| --minstep           | min step in bp [default : 100]                                           |
| --maxdis            | max length of dmr [default : 0]                                          |
| --pvalue            | pvalue cutoff, default: 0.01                                             |
| --FDR               | adjust pvalue cutoff default : 1.0                                       |
| --methdiff          | the cutoff of methylation differention. default: 0.25 [CpG]              |
| --context           | Context for DM. CG/CHG/CHH/C, [C]                                        |
| -h|--help                                                                                      |


DMtools is a BS-Seq methylation file view/process tool based on htslib and libBigWig.

### Single-cell QC (sc-qc)

dmtools supports an ID field in dm files for per-cell tagging. You can summarize basic QC metrics per cell with:

```
dmtools sc-qc -i singlecell.dm -o sc_qc.tsv --context CG --min-coverage 3
```

The output table reports `cell_id`, `n_sites`, `total_coverage`, `mean_coverage`, and coverage-weighted `mean_meth` for each cell ID.

### Single-cell matrix (sc-matrix)

Build a cell × region matrix using the ID field as `cell_id`:

```
dmtools sc-matrix -i singlecell.dm --bed regions.bed -o sc_matrix --context CG --min-coverage 3 --sparse
```

Provide regions via `--bed` (chrom start end [name]) or define fixed windows with `--binsize N`. Sparse Matrix Market output is written by default as `<prefix>.matrix.mtx` along with `barcodes.tsv`, `features.tsv`, and `obs_qc.tsv`; pass `--dense` to write a dense TSV matrix instead.

### Single-cell export (sc-export)

Convert the MatrixMarket bundle into an `.h5ad` file (requires Python with `anndata`, `scipy`, `pandas`):

```
dmtools sc-export -i singlecell.dm --bed regions.bed -o sc_matrix --to h5ad
```

If the Python dependencies are missing, the command still writes the MatrixMarket bundle and prints a warning about installing the missing packages.

### Single-cell aggregation (sc-aggregate)

Aggregate single-cell methylation into group-level profiles using a cell → group mapping:

```
dmtools sc-aggregate -i singlecell.dm --groups cell_to_group.tsv --bed regions.bed -o sc_agg --context CG --min-coverage 3 --dense
```

This writes per-group region summaries and, with `--dense`, a group × region matrix alongside optional feature and group listings.

For more information, please see https://dmtools-docs.rtfd.io/ and the format spec in `docs/spec.md`.

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://dmtools-docs.readthedocs.io/en/latest/function/bam2dm.html

## Debugging and validation

The `bam2dm` converter benefits from running with the address/undefined behaviour sanitizers enabled when investigating platform-specific crashes or corrupted DM output:

```
make clean
CXXFLAGS="-O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer" make bam2dm
```

The `bam2dm` binary built by `make bam2dm` comes from `bam2dm.cpp`; `dmtools bam2dm` shells out to that executable.

Run the converter with verbose diagnostics and on-disk validation enabled to surface silent failures early:

```
./dmtools bam2dm --debug --validate-output \
  -g test/test.chr1.f \
  -b /tmp/minimal.bam \
  -m /tmp/minimal.dm \
  --chunk-by chrom --bin-size 2000 --threads 1
```

Bin mode uses fixed-size bins (default 2000 bp) to distribute work across threads while keeping a single, deterministic writer. Each worker opens its own BAM+index handle and appends results to a per-thread shard under `<methratio>.parts/thread_<tid>.tmp`, capping the number of temporary files to the thread count instead of the bin count. The main thread replays each shard in `(chrom,start)` order into the final dm before running `--validate-output` if requested. A BAM index (`.bai`/`.csi`) is required; missing indexes cause the command to fail fast.

Example bin-mode invocation with validation:

```
./dmtools bam2dm --debug --validate-output \
  -g test/test.chr1.f \
  -b /tmp/minimal.bam \
  -m /tmp/minimal.dm \
  --chunk-by bin --bin-size 2000 --threads 4
```

To compare deterministic output across thread counts, rerun the same input with `--threads 2` (or higher) and ensure the `--validate-output` check succeeds for each run.

To validate an existing DM without rewriting it, use the structural check mode, which returns a non-zero exit code on corruption and prints a short success message otherwise:

```
./dmtools bam2dm --check /tmp/minimal.dm
dmtools validate -i /tmp/minimal.dm
```

### next steps
- [x] add support for NOMe-seq
- [x] add quality contorl and alignment (with fastp/bwa)
- [ ] test suport for scBS
