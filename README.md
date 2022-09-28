# BMtools

Package of C and C++ scripts for generation and procession of DNA Methylation BigWig (bm) files

Compared with bigwig, bm format contains the coverage of DNA methylation sites, methylation context, positive and negative strand and gene ID information. Moreover, bm files still support the visualization of genome browsers as bigwig files.

BMtools provide file format conversion, quick view, methylation level calculation, differential DNA methylation level calculation and other functions for bm format.


| [ bmtools ]         | process with bm file                                                    |
|---------------------|--------------------------------------------------------------------------|
|                     | bmtools <mode> [opnions]                                                 |
| [mode]              | mr2bm view overlap regionstats bodystats profile chromstats             |
| mr2bm              | convert txt meth file to bm format                                      |
| view                | bm format to txt meth                                                   |
| viewheader          | view header of bm file                                                  |
| overlap             | overlap cytosine site with more than two bm files                       |
| regionstats         | calculate DNA methylation level of per region                            |
| bodystats           | calculate DNA methylation level of body, upstream and downstream.        |
| profile             | calculate DNA methylation profile                                        |
| chromstats          | calculate DNA methylation level across chromosome                        |
| -h|--help                                                                                      |

| bmDMR               | detect DMC/DMR based on bm file                                         |
|---------------------|--------------------------------------------------------------------------|
| -p                  | output file prefix                                                       |
| -1                  | sample1 methy bm files, sperate by comma.                               |
| -2                  | sample2 methy bm files, sperate by comma.                               |
| --mindmc            | min dmc sites in dmr region. [default : 4]                               |
| --minstep           | min step in bp [default : 100]                                           |
| --maxdis            | max length of dmr [default : 0]                                          |
| --pvalue            | pvalue cutoff, default: 0.01                                             |
| --FDR               | adjust pvalue cutoff default : 1.0                                       |
| --methdiff          | the cutoff of methylation differention. default: 0.25 [CpG]              |
| --context           | Context for DM. CG/CHG/CHH/C, [C]                                        |
| -h|--help                                                                                      | 


This is a BS-Seq results bm format file view and process tool based on libBigWig.

For more information, please see https://batmeth2-bm.readthedocs.io/en/latest/function/BMtools.html

And calmeth in batmeth2-bm can convert align bs bam file to bm file, https://batmeth2-bm.readthedocs.io/en/latest/function/Calmeth.html

