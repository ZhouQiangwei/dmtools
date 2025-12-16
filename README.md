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

Package of C and C++ scripts for generation and procession of DNA Methylation (DM) files

We present a new and efficient design for storing DNA methylation (DM) data after mapping in compressed binary indexed DM format. Our format significantly reduces storage space by 80%-95% compared to commonly used file formats for DNA methylation data after mapping. To enhance the processing of DNA methylation data in DM format, we have developed DMtools, a comprehensive toolkit that offers utilities such as rapid and random access, computation of DNA methylation profiles across genes, and analysis of differential DNA methylation.

DM format contains the coverage of DNA methylation sites, methylation context, positive and negative strand and cell ID (for single-cell DNA methylation data) information. 

DMtools provide file format conversion, quick view, methylation level calculation, differential DNA methylation level calculation and other functions for dm format.

For more information, please see https://dmtools-docs.rtfd.io

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://dmtools-docs.readthedocs.io/en/latest/function/bam2dm.html


| [ dmtools ]         | process with dm file                                                    |
|---------------------|--------------------------------------------------------------------------|
|                     | dmtools <mode> [opnions]                                                 |
| [mode]              | index align mr2dm view overlap regionstats bodystats profile chromstats |
| index               | build index for genome                                                  |
| align               | alignment fastq                                                         |
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


This is a BS-Seq results dm format file view and process tool based on htslib and libBigWig.

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

Provide regions via `--bed` (chrom start end [name]) or define fixed windows with `--binsize N`. Sparse Matrix Market output is written by default along with `barcodes.tsv` and `features.tsv`; pass `--dense` to write a dense TSV matrix instead.

### Single-cell aggregation (sc-aggregate)

Aggregate single-cell methylation into group-level profiles using a cell → group mapping:

```
dmtools sc-aggregate -i singlecell.dm --groups cell_to_group.tsv --bed regions.bed -o sc_agg --context CG --min-coverage 3 --dense
```

This writes per-group region summaries and, with `--dense`, a group × region matrix alongside optional feature and group listings.

For more information, please see https://dmtools-docs.rtfd.io/

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://dmtools-docs.readthedocs.io/en/latest/function/bam2dm.html

## Debugging and validation

The `bam2dm` converter benefits from running with the address/undefined behaviour sanitizers enabled when investigating platform-specific crashes or corrupted DM output:

```
make clean
CXXFLAGS="-O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer" make bam2dm
```

The `bam2dm` binary built by `make bam2dm` comes from `bam2dm.cpp`; `dmtools bam2dm` shells out to that executable. The alternate source `bam2dm-mp.cpp` is retained for reference but is **not** compiled or invoked by `dmtools`.

Run the converter with verbose diagnostics and on-disk validation enabled to surface silent failures early:

```
./dmtools bam2dm --debug --validate-output \
  -g test/test.chr1.f \
  -b /tmp/minimal.bam \
  -m /tmp/minimal.dm \
  --chunk-by bin --bin-size 2000 --threads 1
```

To compare deterministic output across thread counts, rerun the same input with `--threads 2` (or higher) and ensure the `--validate-output` check succeeds for each run.

To validate an existing DM without rewriting it, use the structural check mode, which returns a non-zero exit code on corruption and prints a short success message otherwise:

```
./dmtools bam2dm --check /tmp/minimal.dm
```

### next steps
- [x] add support for NOMe-seq
- [x] add quality contorl and alignment (with fastp/bwa)
- [ ] test suport for scBS
