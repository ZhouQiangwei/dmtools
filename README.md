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
| [mode]              | mr2dm view overlap regionstats bodystats profile chromstats             |
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

For more information, please see https://dmtools-docs.rtfd.io/ 

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://dmtools-docs.readthedocs.io/en/latest/function/bam2dm.html

### next steps
- [x] add support for NOMe-seq
- [x] add quality contorl and alignment (with fastp/bwa)
- [ ] test suport for scBS
