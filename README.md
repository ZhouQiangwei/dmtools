# DMtools

Package of C and C++ scripts for generation and procession of DNA Methylation (dm) files

Compared with bigwig, dm format contains the coverage of DNA methylation sites, methylation context, positive and negative strand and gene ID information. 

DMtools provide file format conversion, quick view, methylation level calculation, differential DNA methylation level calculation and other functions for dm format.


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

For more information, please see https://batmeth2-dm.readthedocs.io/en/latest/function/DMtools.html

And calmeth in batmeth2-dm can convert align bs bam file to dm file, https://batmeth2-dm.readthedocs.io/en/latest/function/Calmeth.html

