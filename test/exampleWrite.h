#include "bigWig.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

FILE* File_Open(const char* File_Name,const char* Mode);
int main_view_all(bigWigFile_t *fp);
int main_view(bigWigFile_t *fp, char *region);
int main_view_bedfile(char *inbmF, char *bedfile, int type);
int calchromstats(char *outbmfile, char *method, int chromstep, int stepoverlap, uint8_t strand);
int calregionstats(char *outbmfile, char *method, char *region, uint8_t pstrand);
int calregionstats_file(char *outbmfile, char *method, char *bedfile, int format);
int bw_overlap_region(char *inbmF1, char *inbmF2, char *region, uint8_t pstrand);
int bw_overlap_region_mul(char *inbmFs, char *region, uint8_t pstrand);
int bw_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand);
int bw_overlap_file_mul(char *inbmF1, char *bedfile);
int bw_overlap_file(char *inbmF1, char *inbmF2, char *bedfile);
int bw_overlap_all_mul(char *inbmFs, uint8_t pstrand);
int bw_overlap_mul(bigWigFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand);
int bw_overlap(bigWigFile_t *ifp1, bigWigFile_t *ifp2, char *chrom, int start, int end, uint8_t strand);
int calprofile(char *outbmfile, int upstream, int downstream, double profilestep, double profilemovestep, char *bedfile);
int calprofile_gtf(char *outbmfile, int upstream, int downstream, double profilestep, double profilemovestep, char *gtffile, int format);
