# BMtools

Package of C scripts for generation and procession of DNA Methylation BigWig (mbw) files

Compared with bigwig, mbw format contains the coverage of DNA methylation sites, methylation context, positive and negative strand and gene ID information. Moreover, mbw files still support the visualization of genome browsers as bigwig files.

BMtools provide file format conversion, quick view, methylation level calculation, differential DNA methylation level calculation and other functions for mbw format.

+---------------------+--------------------------------------------------------------------------+
| **[ Main paramaters ]**                                                                        |
+=====================+==========================================================================+
|                     | bmtools <mode> [opnions]                                                 |
+---------------------+--------------------------------------------------------------------------+
|Usage:                                                                                          |
+---------------------+--------------------------------------------------------------------------+
| [mode]              | mr2mbw view overlap regionstats bodystats profile chromstats             |
+---------------------+--------------------------------------------------------------------------+
| mr2mbw              | convert txt meth file to mbw format                                      |
+---------------------+--------------------------------------------------------------------------+
| view                | mbw format to txt meth                                                   |
+---------------------+--------------------------------------------------------------------------+
| viewheader          | view header of mbw file                                                  |
+---------------------+--------------------------------------------------------------------------+
| overlap             | overlap cytosine site with more than two mbw files                       |
+---------------------+--------------------------------------------------------------------------+
| regionstats         | calculate DNA methylation level of per region                            |
+---------------------+--------------------------------------------------------------------------+
| bodystats           | calculate DNA methylation level of body, upstream and downstream.        |
+---------------------+--------------------------------------------------------------------------+
| profile             | calculate DNA methylation profile                                        |
+---------------------+--------------------------------------------------------------------------+
| chromstats          | calculate DNA methylation level across chromosome                        |
+---------------------+--------------------------------------------------------------------------+
| -h|--help                                                                                      |
+---------------------+--------------------------------------------------------------------------+

This is a BS-Seq results mbw format file view and process tool.

For more information, please see https://batmeth2-mbw.readthedocs.io/en/latest/function/BMtools.html
