/*
* This tools is used for make bigmeth file and process bigmeth file.
* main function contains:
* bam2bm calculate methleval in loci and region with bam file, output bm file.
* view   view bigmeth as txt format, 20211101
* mr2bm  convert methratio file from batmeth2 to bigmeth, 20211103
* bm2mr  convert bigmeth file to methratio txt file, 20211101
* calregion   calculate DNA methylation level of all regions.
* plotbasic
* plotprofile
* plotheatmap
* plotchrom
* SNV   SNP and INDEL detection.
* realign   realign of bsalign bam file.
* 
* make clean && make && gcc libBigWig.a libBigWig.so test/exampleWrite.c -o exampleWrite
* ./exampleWrite mr2bm -g ~/practice/Genome/hg38/hg38.chr.fa.len -E -C -m test.f -S --Cx -o test.bm -r chr1:0-100,chr1:16766-16890
* gcc test/exampleWrite.c -o exampleWrite -I. -L. -lBigWig -Wl,-rpath /public/home/qwzhou/software_devp/batmeth2-bma/src/bmtools/
*/
#include "binaMeth.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <pthread.h>

FILE* File_Open(const char* File_Name,const char* Mode);
char *strand_str[] = {"+", "-", "."};
char *context_str[] = {"C", "CG", "CHG", "CHH"};
int main_view_all(binaMethFile_t *fp, FILE* outfileF, char *outformat, binaMethFile_t *ofp, char* filterchrom);
int main_view(binaMethFile_t *ifp, char *region, FILE* outfileF, char *outformat, binaMethFile_t *ofp);
int main_view_bedfile(char *inbmF, char *bedfile, int type, FILE* outfileF, char *outformat, binaMethFile_t *ofp);
int calchromstats(char *inbmfile, char *method, int chromstep, int stepoverlap, uint8_t strand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calregionstats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calregionstats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calbodystats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calbodystats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int bm_overlap_region(char *inbmF1, char *inbmF2, char *region, uint8_t pstrand);
int bm_overlap_region_mul(char *inbmFs, char *region, uint8_t pstrand);
int bm_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand);
int bm_overlap_file_mul(char *inbmF1, char *bedfile);
int bm_overlap_file(char *inbmF1, char *inbmF2, char *bedfile);
int bm_overlap_all_mul(char *inbmFs, uint8_t pstrand);
int bm_overlap_mul(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand);
int bm_overlap(binaMethFile_t *ifp1, binaMethFile_t *ifp2, char *chrom, int start, int end, uint8_t strand);
int calprofile(char *inbmfile, int upstream, int downstream, double profilestep, double bodyprofilestep, double bodyprofilemovestep, double profilemovestep, char *bedfile, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX, uint8_t filt_strand);
int calprofile_gtf(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *gtffile, int format, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX, uint8_t filt_strand);
void getcontext(uint8_t context, char* str_context);
void delete_char2(char* str,char target,char target2);
void calbody_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF);
void calregion_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* strand, FILE* outfileF);
void calregion_weighted_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, uint32_t *countC, uint32_t *countCT);
int main_view_file(binaMethFile_t *ifp, char *bedfile, FILE* outfileF, char *outformat, binaMethFile_t *ofp);
void bmfileinit(binaMethFile_t *ofp, binaMethFile_t *ifp, char* outfile, int zoomlevel);
void bmPrintHdr(binaMethFile_t *bm);
void bmPrintIndexNode(bmRTreeNode_t *node, int level);
char *fastStrcat(char *s, char *t);
void onlyexecuteCMD(const char *cmd, const char *errorinfor);
size_t get_executable_path( char* processdir,char* processname, size_t len);
unsigned long get_chr_len(binaMethFile_t *bm, char* chrom);
int main_view_bm(binaMethFile_t *ifp, char *region, FILE* outfileF, char *outformat, binaMethFile_t *ofp, \
    char** chromsUse, uint32_t* starts, uint32_t* ends, float* values, uint16_t* coverages, uint8_t* strands, \
    uint8_t* contexts, char** entryid);

#define MAX_LINE_PRINT 1000100
#define MAX_BUFF_PRINT 20000000
const char* Help_String_main="Command Format :  dmtools <mode> [opnions]\n"
		"\nUsage:\n"
        "\t  [mode]         bam2dm mr2dm view viewheader overlap regionstats bodystats profile chromstats\n\n"
        "\t  bam2dm         calculate DNA methylation (BM format) with BAM file\n"
        "\t  mr2dm          convert txt meth file to dm format\n"
        "\t  view           dm format to txt/dm meth\n"
        "\t  viewheader     view header of dm file\n"
        "\t  overlap        overlap cytosine site with more than two dm files\n"
        "\t  regionstats    calculate DNA methylation level of per region\n"
        "\t  bodystats      calculate DNA methylation level of body, upstream and downstream.\n"
        "\t  profile        calculate DNA methylation profile\n"
        "\t  chromstats     calculate DNA methylation level across chromosome\n"
        "\t  addzm          add or change zoom levels for dm format, need for browser visulization\n"
        "\t  stats          coverage and methylation level distribution of data\n";

const char* Help_String_bam2dm="Command Format :  dmtools bam2dm [options] -g genome.fa -i/-b <SamfileSorted/BamfileSorted> -m <methratio dm outfile prefix>\n"
		"\nUsage: dmtools bam2dm -g genome.fa -b align.sort.bam -m meth.dm\n"
        "\t [bam2dm] mode paramaters, required\n"
		"\t-g|--genome           genome fasta file\n"
		"\t-b|--binput           Bam format file, sorted by chrom.\n"
        "\t-m|--methratio        Prefix of methratio.dm output file\n"
		"\t [bam2dm] mode paramaters, options\n"
        "\t-n|--Nmismatch        Number of mismatches, default 0.06 percentage of read length. [0-1]\n"
		"\t-Q                    caculate the methratio while read QulityScore >= Q. default:20\n"
		"\t-c|--coverage         >= <INT> coverage. default:4\n"
		"\t--maxcoverage         <= <INT> coverage. default:500\n"
		"\t-nC                   >= <INT> nCs per region. default:1\n"
		"\t-r|--remove_dup       REMOVE_DUP, default:false\n"
        "\t--mrtxt               print prefix.methratio.txt file\n"
        "\t--zl                  The maximum number of zoom levels. [0-10], default: 2\n"
        "\t-i|--input            Sam format file, sorted by chrom.\n"
        "\t-C                    print coverage in DM file\n"
        "\t-S                    print strand in DM file\n"
        "\t--Cx                  print context in DM file\n"
        "\t-E                    print end in DM file\n"
        "\t--Id                  print ID in DM file\n"
        "\t-h|--help";

const char* Help_String_mr2dm="Command Format :  dmtools mr2dm [opnions] -g genome.fa.fai -m methratio.txt -o outmeth.dm\n"
		"\nUsage: dmtools mr2dm -C -S --Cx -g genome.fa.fai -m meth.txt -o meth.dm\n"
        "\t [mr2dm] mode paramaters, required\n"
		"\t-g                    chromosome size file.\n"
        "\t-m                    methratio file\n"
        "\t-o|--outdm            output DM file\n"
        "\t [mr2dm] mode paramaters, options\n"
        "\t-C                    print coverage\n"
        "\t-S                    print strand\n"
        "\t--Cx                  print context\n"
        "\t-E                    print end\n"
        "\t--Id                  print ID\n"
        "\t--CF                  coverage filter, >=[int], default 4.\n"
        "\t--sort Y/N            make chromsize file and meth file in same coordinate, default Y\n"
        "\t--zl                  The maximum number of zoom levels. [0-10]\n"
        "\t-f                    file format. methratio, bedmethyl, bismark or bedsimple [default methratio]\n"
        "\t  methratio           chrom start strand context meth_reads cover\n"
        "\t  bedmethyl           chrom start end name * strand * * * coverage meth_reads\n"
        "\t  bismark             chrom start strand coverC coverT context\n"
        "\t  bedsimple           chrom start end id strand context meth_reads coverage\n"
        "\t--pcontext            CG/CHG/CHH/C, needed when bedmethyl format, default C\n"
//        "\t--context             [0/1/2/3] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context.\n"
        "\t--fcontext            CG/CHG/CHH/ALL, only convert provide context in methratio file or bedsimple, default ALL\n"
        "\t--sv                  add strand information to meth level. eg 0.5, -0.5\n"
        "\tNote. meth ratio file must be sorted by chrom and coordinate. ex. sort -k1,1 -k2,2n\n"
		"\t-h|--help";

const char* Help_String_view="Command Format :  dmtools view [opnions] -i meth.dm\n"
		"\nUsage:\n"
        "\t [view] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t [view] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by space. chr1:1-2900 chr2:1-200\n"
        "\t--chr                 chromosome for view\n"
        "\t--bed                 bed file for view, format: chrom start end (strand).\n"
        "\t--strand              [0/1/2] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information\n"
        "\t--context             [0/1/2/3] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context.\n"
        "\t--mincover            >= minumum coverage show, default: 0\n"
        "\t--maxcover            <= maximum coverage show, default: 10000\n"
        "\t--outformat           txt or dm format [txt]\n"
        "\t--zl                  The maximum number of zoom levels. [0-10], valid for dm out\n"
		"\t-h|--help";

const char* Help_String_stats="Command Format :  dmtools stats [opnions] -i meth.dm\n"
        "\nUsage:\n"
        "\t [stats] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t [stats] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t--tc                  total number of cytosine and guanine in genome, we will use the number of site in dm file if you not provide --tc\n"
        "\t-r                    region for calculate stats, can be seperated by space. chr1:1-2900 chr2:1-200\n"
        "\t--bed                 bed file for calculate stats, format: chrom start end (strand).\n"
        "\t--strand              [0/1/2] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information\n"
        "\t--context             [0/1/2/3] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context.\n"
        "\t-s                    size of bin for count cytosine number.\n"
        "\t--mincover            >= minumum coverage show, default: 0\n"
        "\t--maxcover            <= maximum coverage show, default: 10000\n"
        "\t-h|--help";

const char* Help_String_addzm="Command Format :  dmtools addzm [opnions] -i meth.dm -o meth.zm2.dm\n"
		"\nUsage: add zoom levels for dm\n"
        "\t [addzm] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t [addzm] mode paramaters, options\n"
        "\t-o                    output dm file\n"
        "\t--zl                  The maximum number of zoom levels. [0-10], valid for dm out\n"
        "\t-r                    region for view, can be seperated by space. chr1:1-2900 chr2:1-200;\n"
        "\t--bed                 bed file for view, format: chrom start end (strand).\n"
        "\t--strand              [0/1/2] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information\n"
        "\t--context             [0/1/2/3] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context.\n"
        "\t--mincover            >= minumum coverage show, default: 0\n"
        "\t--maxcover            <= maximum coverage show, default: 10000\n"
		"\t-h|--help";

const char* Help_String_viewheader="Command Format :  dmtools viewheader -i meth.dm\n"
		"\nUsage:\n"
        "\t [view] mode paramaters, required\n"
        "\t-i                    input DM file\n"
		"\t-h|--help";

const char* Help_String_overlap="Command Format :  dmtools overlap [opnions] -i meth1.dm -i2 meth2.dm\n"
		"\nUsage:\n"
        "\t [overlap] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t-i2                   input DM file2\n"
        "\t [overlap] mode paramaters, options\n"
        //"\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by space. chr1:1-2900 chr2:1-200\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t [overlap] mode paramaters\n"
        "\t--dmfiles             input DM files, seperated by comma. This is no need if you provide -i and -i2.\n"
		"\t-h|--help";

const char* Help_String_regionstats="Command Format :  dmtools regionstats [opnions] -i dm --bed bedfile\n"
		"\nUsage:\n"
        "\t [regionstats] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
        "\t [regionstats] mode paramaters, options\n"
        "\t-o                    output prefix [stdout]\n"
        "\t-r                    region for view, can be seperated by space. chr1:1-2900 chr2:1-200,+\n"
        "\t--method              weighted/ mean\n"
        "\t--strand              [0/1/2/3] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information, 3 calculate and print strand meth level seperately\n"
        "\t--context             [0/1/2/3/4] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context, 4 calculate and print strand meth level seperately\n"
        "\t--printcoverage       [0/1] print countC and coverage instead of methratio. [0]\n"
        "\t--print2one           [int] print all the countC and coverage results of C/CG/CHG/CHH context methylation to same file, only valid when --printcoverage 1. 0 for no, 1 for yes. [0]\n"
//        "\t--mincover            >= minumum coverage show, default: 0\n"
//        "\t--maxcover            <= maximum coverage show, default: 10000\n"
		"\t-h|--help";

const char* Help_String_bodystats="Command Format :  dmtools bodystats [opnions] -i dm --bed bedfile\n"
		"\nUsage:\n"
        "\t [bodystats] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
        "\t [bodystats] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by space. chr1:1-2900 chr2:1-200,+\n"
        "\t--method              weighted/ mean\n"
        "\t--regionextend        also calculate DNA methylation level of upstream and downstream N-bp window. default 2000.\n"
        "\t--strand              [0/1/2/3] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information, 3 calculate and print strand meth level seperately, only valid while -r para\n"
        "\t--context             [0/1/2/3/4] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context, 4 calculate and print strand meth level seperately, default: 4.\n"
		"\t--printcoverage       [0/1] also print countC and coverage of body instead of methratio. [0]\n"
        "\t--print2one           [int] print all the countC and coverage results of C/CG/CHG/CHH context methylation to same file, only valid when --printcoverage 1. 0 for no, 1 for yes. [0]\n"
        "\t-h|--help";

const char* Help_String_profile="Command Format :  dmtools profile [opnions] -i meth.dm --bed bedfile\n"
		"\nUsage:\n"
        "\t [profile] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
		"\t [profile] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-p                    [int] threads\n"
        "\t--regionextend        region extend for upstream and downstram, [2000]\n"
        "\t--profilestep         [double] step mean bin size for chromosome region, default: 0.02 (2%)\n"
        "\t--profilemovestep     [double] step move, default: profilestep/2, if no overlap, please define same as --profilestep\n"
//        "\t--bodyprofilestep     [double] step mean bin size for element body region, default: 0.02 (2%)\n"
//        "\t--bodyprofilemovestep [double] step move for body region, default: 0.01, if no overlap, please define same as --bodyprofilestep\n"
        "\t--profilemode         calculate the methylation matrix mode of every region or gene. 0 for gene and flanks mode, 1 for tss and flanks, 2 for tts and flanks, 3 for region center and flanks. [0]\n"
        "\t--bodyX               [double] the gene body bin size is N multiple of the bin size of the upstream and downstream extension region. [1]\n"
        "\t--matrixX             [int] the bin size is N times of the profile bin size, so the bin size should be N*profilestep [5], please note N*profilestep must < 1 and N must >= 1, used for calculation per gene.\n"
        "\t--strand              [0/1/2] strand for calculate, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information, [2]\n"
        "\t--context             [0/1/2/3] context for calculate, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context. [0]\n"
        "\t--print2one           [int] print all the matrix results of C/CG/CHG/CHH context methylation to same file. 0 for no, 1 for yes. [0]\n"
        "\t-h|--help";

const char* Help_String_chromstats="Command Format :  dmtools chromstats [opnions] -i meth.dm\n"
		"\nUsage:\n"
        "\t [chromstats] mode paramaters, required\n"
        "\t-i                    input DM file\n"
        "\t [chromstats] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t--method              weighted/ mean/ min/ max/ cover/ dev\n"
        "\t--chromstep           [int] step mean bin size for chromosome region, default: 100000\n"
        "\t--stepmove            [int] step move, default: 50000, if no overlap, please define same as --chromstep\n"
        "\t--context             [0/1/2/3/4] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context, 4 calculate and print context meth level seperately. [4]\n"
        "\t--fstrand             [0/1/2/3] strand for calculation, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' without strand information, 3 calculate and print strand meth level seperately. [2]\n"
		"\t--printcoverage       [0/1] print countC and coverage instead of methratio. [0]\n"
        "\t--print2one           [int] print all the countC and coverage results of C/CG/CHG/CHH context methylation to same file, only valid when --printcoverage 1. 0 for no, 1 for yes. [0]\n"
        "\t-h|--help";

int regionextend = 2000;
//0, c; 1 cg; 2 chg; 3 chh; 4 seperate
uint8_t filter_context = 4;
//0 +; 1 -; 2 . all;
uint8_t filter_strand = 2;
int mincover = 0;
int maxcover = 10000;
int printcoverage = 0; // print countC and countCT instead of methratio
int NTHREAD = 10;
int *Fcover;
unsigned long *mPs;
//mP1 = 0, mP2 = 0, mP3 = 0, mP4 = 0, mP5 = 0;
int statsSize = 10;
int alwaysprint = 0;

int main(int argc, char *argv[]) {
    binaMethFile_t *fp = NULL;
    char *filterchrom = NULL;
    char **chroms = (char**)malloc(sizeof(char*)*MAX_LINE_PRINT);
    if(!chroms) goto error;
    char **chromsUse = malloc(sizeof(char*)*MAX_LINE_PRINT);
    char **entryid = malloc(sizeof(char*)*MAX_LINE_PRINT);
    uint32_t *chrLens = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    uint32_t *starts = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    uint32_t *ends = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    float *values = malloc(sizeof(float) * MAX_LINE_PRINT);
    uint16_t *coverages = malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
    uint8_t *strands = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
    uint8_t *contexts = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);

    char* chromlenf = malloc(100*sizeof(char));
    int chromlenf_yes = 0;
    char* outformat = malloc(100); strcpy(outformat, "txt");
    Fcover = malloc(sizeof(int)*16);
    unsigned long totalCG = 0;
    int i = 0;
    uint32_t write_type = 0x8000;
    char methfile[100]; char *outbmfile = NULL;
    char *inbmfile = malloc(100);
    char *bmfile2 = malloc(100);
    char mode[10]; char *region = NULL;
    char *inbmfiles = malloc(300);
    int inbm_mul = 0;
    char *method = malloc(100);
    strcpy(method, "weighted");
    int chromstep = 100000;
    int stepoverlap = 50000;
    
    char *bedfile = NULL;
    char *gtffile = NULL;
    char *gfffile = NULL;
    int strandmeth = 0;

    char *mrformat = malloc(100);
    strcpy(mrformat, "methratio");
    char *pcontext = malloc(10);
    strcpy(pcontext, "C");
    char *filtercontext = malloc(10);
    strcpy(filtercontext, "ALL");
    // for gene file
    int upstream = 2000, downstream = 2000;
    double profilestep = 0.02, profilemovestep = 0.01; int defineprofilemovestep = 0;
    double bodyprofilestep = 0.02, bodyprofilemovestep = 0.01;
    unsigned int mcover_cutoff = 4; unsigned long TotalC = 0;
    int zoomlevel = 2;
    char *sortY = malloc(10); strcpy(sortY, "Y");
    char* outfile = NULL;
    int profilemode = 0; //0 gene and flanks, 1 center and flanks, 2 tss and flanks
    char *profilemode_str = malloc(sizeof(char)*10);
    double bodyX = 1; // N times of bin size for gene body
    int matrixX = 5; int print2one = 0; // print all matrix to one file.
    if(argc>3){
        strcpy(mode, argv[1]);
        // mr2bam view overlap region
    }else if(argc>1){
        strcpy(mode, argv[1]);
        //mr2bm/ view/ overlap/ regionstats/ profile/ chromstats
        if(strcmp(mode, "mr2dm") == 0){
           fprintf(stderr, "%s\n", Help_String_mr2dm); 
        }else if(strcmp(mode, "view") == 0){
           fprintf(stderr, "%s\n", Help_String_view); 
        }else if(strcmp(mode, "stats") == 0){
           fprintf(stderr, "%s\n", Help_String_stats);
        }else if(strcmp(mode, "viewheader") == 0){
           fprintf(stderr, "%s\n", Help_String_viewheader); 
        }else if(strcmp(mode, "overlap") == 0){
           fprintf(stderr, "%s\n", Help_String_overlap); 
        }else if(strcmp(mode, "regionstats") == 0){
           fprintf(stderr, "%s\n", Help_String_regionstats); 
        }else if(strcmp(mode, "bodystats") == 0){
           fprintf(stderr, "%s\n", Help_String_bodystats); 
        }else if(strcmp(mode, "profile") == 0){
           fprintf(stderr, "%s\n", Help_String_profile); 
        }else if(strcmp(mode, "chromstats") == 0){
           fprintf(stderr, "%s\n", Help_String_chromstats); 
        }else if(strcmp(mode, "bam2dm") == 0){
           fprintf(stderr, "%s\n", Help_String_bam2dm); 
        }else if(strcmp(mode, "addzm") == 0){
           fprintf(stderr, "%s\n", Help_String_addzm); 
        }else{
            fprintf(stderr, "Please define correct mode!!!\n");
            fprintf(stderr, "%s\n", Help_String_main);
        }
        exit(0);
    }else{
        //fprintf(stderr, "Please define mode!!!\n");
        fprintf(stderr, "%s\n", Help_String_main);
        exit(0);
    }
    char bam2bm_paras[1024];
    if(strcmp(mode, "bam2dm") == 0){
        for(i=1; i< argc; i++){
            strcat(bam2bm_paras, argv[i]);
            strcat(bam2bm_paras, " ");
        }
    }else{
        for(i=1; i< argc; i++){
            if(strcmp(argv[i], "-g") == 0){
                strcpy(chromlenf, argv[++i]);
                chromlenf_yes++;
            }else if(strcmp(argv[i], "-C") == 0){
                write_type |= BM_COVER;
            }else if(strcmp(argv[i], "-S") == 0){
                write_type |= BM_STRAND;
            }else if(strcmp(argv[i], "--Cx") == 0){
                write_type |= BM_CONTEXT;
            }else if(strcmp(argv[i], "--Id") == 0){
                write_type |= BM_ID;
            }else if(strcmp(argv[i], "-E") == 0){
                write_type |= BM_END;
            }else if(strcmp(argv[i], "-m") == 0){
                strcpy(methfile, argv[++i]);
            }else if(strcmp(argv[i], "--outdm") == 0){
                outbmfile = malloc(100);
                strcpy(outbmfile, argv[++i]);
            }else if(strcmp(argv[i], "--outformat") == 0){
                strcpy(outformat, argv[++i]);
            }else if(strcmp(argv[i], "-i") == 0){
                strcpy(inbmfile, argv[++i]);
            }else if(strcmp(argv[i], "--CF") == 0){
                mcover_cutoff = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--dmfiles") == 0){
                strcpy(inbmfiles, argv[++i]);
                inbm_mul = 1;
            }else if(strcmp(argv[i], "-i2") == 0){
                strcpy(bmfile2, argv[++i]);
            }else if(strcmp(argv[i], "-r") == 0){
                region = malloc(1000);
                strcpy(region, argv[++i]);
                char* rtemp = malloc(200);
                while(i+1<argc && argv[++i][0]!='-'){
                    strcpy(rtemp, argv[i]);
                    strcat(region, ";");
                    strcat(region, rtemp);
                }
                i--;
                free(rtemp);
            }else if(strcmp(argv[i], "--method") == 0){
                strcpy(method, argv[++i]);
            }else if(strcmp(argv[i], "--chromstep") == 0){
                chromstep = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--chrom") == 0 || strcmp(argv[i], "--chr") == 0){
                filterchrom = malloc(200);
                strcpy(filterchrom, argv[++i]);
            }else if(strcmp(argv[i], "--stepmove") == 0){
                stepoverlap = atoi(argv[++i]);
            }else if(strcmp(argv[i], "-p") == 0){
                NTHREAD = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--tc") == 0){
                totalCG = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--profilestep") == 0){
                profilestep = atof(argv[++i]);
                assert(profilestep>0 && profilestep<1);
            }else if(strcmp(argv[i], "--regionextend") == 0){
                regionextend = atoi(argv[++i]);
                assert(regionextend>0);
            }else if(strcmp(argv[i], "--profilemovestep") == 0){
                defineprofilemovestep = 1;
                profilemovestep = atof(argv[++i]);
            }else if(strcmp(argv[i], "--bodyprofilestep") == 0){
                bodyprofilestep = atof(argv[++i]);
                assert(bodyprofilestep>0 && bodyprofilestep<1);
            }else if(strcmp(argv[i], "--bodyprofilemovestep") == 0){
                bodyprofilemovestep = atof(argv[++i]);
            }else if(strcmp(argv[i], "--profilemode") == 0){
                profilemode = atoi(argv[++i]);
                assert(profilemode>=0 && profilemode<4);
                if(profilemode == 0){
                    strcpy(profilemode_str, ".across");
                }else if(profilemode == 1){
                    strcpy(profilemode_str, ".tss");
                }else if(profilemode == 2){
                    strcpy(profilemode_str, ".tts");
                }else if(profilemode == 3){
                    strcpy(profilemode_str, ".center");
                }
            }else if(strcmp(argv[i], "--print2one") == 0){
                print2one = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--printcoverage") == 0){
                printcoverage = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--bodyX") == 0){
                bodyX = atof(argv[++i]);
            }else if(strcmp(argv[i], "--matrixX") == 0){
                matrixX = atoi(argv[++i]);
                assert(matrixX>=1);
            }else if(strcmp(argv[i], "--fstrand") == 0 || strcmp(argv[i], "--strand") == 0){
                filter_strand = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--bed") == 0){
                bedfile = malloc(100);
                strcpy(bedfile, argv[++i]);
            }else if(strcmp(argv[i], "--gtf") == 0){
                gtffile = malloc(100);
                strcpy(gtffile, argv[++i]);
            }else if(strcmp(argv[i], "--gff") == 0){
                gfffile = malloc(100);
                strcpy(gfffile, argv[++i]);
            }else if(strcmp(argv[i], "-f") == 0){
                strcpy(mrformat, argv[++i]);
            }else if(strcmp(argv[i], "--pcontext") == 0){
                strcpy(pcontext, argv[++i]);
            }else if(strcmp(argv[i], "--context") == 0){
                filter_context = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--sort") == 0){
                strcpy(sortY, argv[++i]);
            }else if(strcmp(argv[i], "-s") == 0){
               statsSize = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--ap") == 0){
               alwaysprint = 1;
            }else if(strcmp(argv[i], "--zl") == 0){
                zoomlevel = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--fcontext") == 0){
                strcpy(filtercontext, argv[++i]);
            }else if(strcmp(argv[i], "--mincover") == 0){
                mincover = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--maxcover") == 0){
                maxcover = atoi(argv[++i]);
            }else if(strcmp(argv[i], "--sv") == 0){
                strandmeth = 1;
            }else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--out") == 0){
                if(strcmp(mode, "mr2dm") == 0){
                    if(!outbmfile){
                        outbmfile = malloc(100);
                        strcpy(outbmfile, argv[++i]);
                    }
                }else{
                    outfile = malloc(sizeof(char)*100);
                    strcpy(outfile, argv[++i]);
                }
            }
        }
        if(defineprofilemovestep == 0) profilemovestep = profilestep/2;
        upstream = regionextend; downstream = regionextend;
    }

    mPs = calloc(statsSize, sizeof(unsigned long));
    if(strcmp(mode, "bam2dm") == 0){
        fprintf(stderr, "calculate DNA methylation level with dm format\n");
        //exe location
        char processname[1024];
        char abspathtmp[1024];
        get_executable_path(abspathtmp, processname, sizeof(abspathtmp));
        char cmd[1024];
        strcat(cmd, "bam2dm ");
        strcat(cmd, bam2bm_paras);
        onlyexecuteCMD(cmd, Help_String_bam2dm);
        return;
    }
    else if(strcmp(mode, "mr2dm") == 0){
        fprintf(stderr, "mr file format %s\n", mrformat);
        if(chromlenf_yes==0){
            fprintf(stderr, "please provide chrome size file with -g paramater\n");
            exit(0);
        }
        FILE *methF = File_Open(methfile, "r"); 
        char *chrom = malloc(50); char *old_chrom = malloc(50);
        //int MAX_CHROM = 10000;
        char *PerLine = malloc(2000); 
        //char **chromsArray = malloc(sizeof(char*)*MAX_CHROM);
        unsigned long chrprintL = 0;
        if(strcmp(sortY, "Y") == 0){
            fprintf(stderr, "obtained chromosome order in meth ratio file ... \n");
            while(fgets(PerLine,2000,methF)!=0){
                if(PerLine[0] == '#') continue; // remove header #
                //fprintf(stderr, "%s\n", PerLine);
                if(strcmp(mrformat, "methratio") == 0 || strcmp(mrformat, "bedmethyl") == 0 || strcmp(mrformat, "bedsimple") == 0 || strcmp(mrformat, "bismark") == 0){
                    sscanf(PerLine, "%s", chrom);
                }else{
                    fprintf(stderr, "Unexpected mr file format!!!\n");
                    exit(0);
                }

                if(strcmp(old_chrom, chrom)!=0){
                    fprintf(stderr, "chromx %s\n", chrom);
                    chroms[chrprintL++] = strdup(chrom);
                    strcpy(old_chrom, chrom);
                }
            }
        }
        fclose(methF);
        // open chrom size
        FILE* ChromF=File_Open(chromlenf,"r");
        int chrlen = 0, printL = 0, i = 0, inmr = 0, printedchr = 0;
        while(fgets(PerLine,2000,ChromF)!=NULL){
            if(PerLine[0] == '#') continue;
            sscanf(PerLine, "%s%d", chrom, &chrlen);
            //chrLens[printL]=chrlen;
            //chroms[printL++] = strdup(chrom);
            inmr = 0;
            for(i=0;i<chrprintL;i++){
                if(strcmp(chroms[i], chrom) == 0){
                    fprintf(stderr, "chrom %s %d\n", chroms[i], i);
                    chrLens[i] = chrlen;
                    inmr = 1;
                }
            }
            if(inmr == 0){
                chrLens[chrprintL+printedchr]=chrlen;
                chroms[chrprintL+printedchr] = strdup(chrom);
                printedchr++;
            }
            printL++;
        }
        fclose(ChromF);

        fprintf(stderr, "[Mode] ------------- %s\n", mode);
        if(bmInit(1<<17) != 0) {
            fprintf(stderr, "Received an error in bmInit\n");
            return 1;
        }

        fp = bmOpen(outbmfile, NULL, "w");
        fp->type = write_type;
        if(!fp) {
            fprintf(stderr, "An error occurred while opening example_output.dm for writingn\n");
            return 1;
        }
        
        //Allow up to 10 zoom levels, though fewer will be used in practice
        if(bmCreateHdr(fp, zoomlevel)) {
            fprintf(stderr, "== bmCreateHdr ==\n");
            goto error;
        }

        //Create the chromosome lists
        fp->cl = bmCreateChromList(chroms, chrLens, printL); //2
        if(!fp->cl) {
            fprintf(stderr, "== bmCreateChromList ==\n");
            goto error;
        }

        //Write the header
        if(bmWriteHdr(fp)) {
            fprintf(stderr, "== bmWriteHdr ==\n");
            goto error;
        }
        
        //Some example methlevel
        if(DEBUG>1) fprintf(stderr, "====HHH type %d\n", fp->type);
        methF = File_Open(methfile, "r");
        strcpy(old_chrom, "NN"); 
        char *strand = malloc(2), *context = malloc(10), *nameid = malloc(100);
        unsigned start=0, end = 0; unsigned int coverC, coverT, coverage=0; float value=0;
        printL = 0;
        strcpy(context, pcontext);
        char *decide = malloc(10);
        while(fgets(PerLine,2000,methF)!=0){
            if(PerLine[0] == '#') continue; // remove header #
            //fprintf(stderr, "%s\n", PerLine);
            if(strcmp(mrformat, "methratio") == 0){
                sscanf(PerLine, "%s%u%s%s%u%u", chrom, &start, strand, context, &coverC, &coverage);
                //fprintf(stderr, "%f\n", ((double) coverC/ coverage));
                //sprintf(decide, "%.5f", ((double) coverC/ coverage) );
                //value = atof(decide);
                value = ((double) coverC/ coverage);
                //fprintf(stderr, "%f\n", value);
            }else if(strcmp(mrformat, "bedmethyl") == 0){
                //bedmethyl2bm
                //chrom start end name score strand thickStart thickEnd itemRgb reads_coverage percent
                //chr11   61452   61453   .       14      +       61452   61453   0,255,0 14      0
                sscanf(PerLine, "%s%u%u%s%*s%s%*s%*s%*s%u%f", chrom, &start, &end, nameid, strand, &coverage, &value);
                value = value/100;
                coverC = (int)(value*coverage+0.5);
                // cause start from 0
                start++; end++;
            }else if(strcmp(mrformat, "bedsimple") == 0){
                //chrom start end id strand context meth_reads coverage
                sscanf(PerLine, "%s%u%u%s%s%s%u%u", chrom, &start, &end, nameid, strand, context, &coverC, &coverage);
                sprintf(decide, "%0.001f", ((double) coverC/ coverage) );
                value = atof(decide);
            }else if(strcmp(mrformat, "bismark") == 0){
                //<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
                sscanf(PerLine, "%s%u%s%u%u%s", chrom, &start, strand, &coverC, &coverT, context);
                coverage = coverC + coverT;
                value = ((double) coverC/ coverage);
            }else{
                fprintf(stderr, "Unexpected mr file format!!!\n");
                exit(0);
            }
            if(filterchrom && strcmp(filterchrom, chrom) !=0 ) {
                if(TotalC>0) break;
                continue;
            }
            if(coverage<mcover_cutoff) continue;
            if(strcmp(filtercontext, "ALL") != 0 && strcmp(filtercontext, context) != 0){
                continue;
            }
            if(strandmeth == 1) {
                if(strand[0] == '-') value = -value;
            }

            if(strcmp(old_chrom, chrom)!=0){
                if(printL > 1){
                    if(bmAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1))  {
                        fprintf(stderr, "bmAppendIntervals -1\n");
                        goto error;
                    }
                    fprintf(stderr, "Processed %d cytosine sites.\n", printL);
                }
                
                TotalC+=printL;
                printL = 0;
                // init
                chromsUse[printL] = strdup(chrom);
                starts[printL] = start;
                if(end == 0){
                    ends[printL] = start+1;
                }else{
                    ends[printL] = end;
                }
                coverages[printL] = coverage;
                //if(value!=0) {
                //    int tmpv = (int)(value*1000);
                //    value = ((double)tmpv/1000);
                //}
                values[printL] = value;
                
                switch(strand[0]){
                    case '+':
                        strands[printL] = 0;
                        break;
                    case '-':
                        strands[printL] = 1;
                        break;
                    default:
                        strands[printL] = 2;
                        break;
                }
                if(strcmp(context, "C") == 0){
                    contexts[printL] = 0;
                }else if(strcmp(context, "CG") == 0){
                    contexts[printL] = 1;
                }else if(strcmp(context, "CHG") == 0){
                    contexts[printL] = 2;
                }else if(strcmp(context, "CHH") == 0){
                    contexts[printL] = 3;
                }
                if(DEBUG>-1) fprintf(stderr,"## %d start %s %d %d %f %d %d %d\n", printL, chromsUse[printL], starts[printL], ends[printL], values[printL], coverages[printL], strands[printL], context[printL]);
                int response = bmAddIntervals(fp, chromsUse, starts, ends, values, coverages, strands, contexts, 
                entryid, 1);
                fprintf(stderr, "Processing %s chromosome.\n", chrom);
                if(response) {
                    fprintf(stderr, "bmAddIntervals 0\n");
                    goto error;
                }
                strcpy(old_chrom, chrom);
                printL++;
            }else{
                //sscanf(PerLine, "%s%d%s%s%d%d%f", chrom, &starts[printL], strand, context, &coverC, &coverages[printL], &values[printL]);
                //chromsUse[printL] = strdup(chrom);
                //ends[printL] = starts[printL]+1;
                //fprintf(stderr, "\n--222--- %s %d %d %s %s %f\n", chrom, start, coverage, context, strand, value);
                chromsUse[printL] = strdup(chrom);
                starts[printL] = start;
                //ends[printL] = start+1;
                if(end == 0){
                    ends[printL] = start+1;
                }else{
                    ends[printL] = end;
                }
                coverages[printL] = coverage;
                //if(value!=0) {
                //    int tmpv = (int)(value*1000);
                //    value = ((double)tmpv/1000);
                //}
                values[printL] = value;

                switch(strand[0]){
                    case '+':
                        strands[printL] = 0;
                        break;
                    case '-':
                        strands[printL] = 1;
                        break;
                    default:
                        strands[printL] = 2;
                        break;
                }
                if(strcmp(context, "C") == 0){
                    contexts[printL] = 0;
                }else if(strcmp(context, "CG") == 0){
                    contexts[printL] = 1;
                }else if(strcmp(context, "CHG") == 0){
                    contexts[printL] = 2;
                }else if(strcmp(context, "CHH") == 0){
                    contexts[printL] = 3;
                }
                strcpy(old_chrom, chrom);
                printL++;
            }
            if(printL>MAX_LINE_PRINT){
                //We can continue appending similarly formatted entries
                //N.B. you can't append a different chromosome (those always go into different
                if(bmAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) {
                    fprintf(stderr, "bmAppendIntervals 2\n");
                    goto error;
                }
                TotalC+=printL;
                printL = 0;
            }
        } // end read me file
        if(printL > 1) {
            if(DEBUG>1) fprintf(stderr, "--last print %d %d %d\n", starts[printL-1], ends[printL-1], printL-1);
            if(bmAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) {
                fprintf(stderr, "bmAppendIntervals 3\n");
                goto error;
            }
            TotalC+=printL;
            printL = 0; 
        }
        fprintf(stderr, "\nprocess %ld cytosine site\n", TotalC);
        //Add a new block of entries with a span. Since bmAdd/AppendIntervals was just used we MUST create a new block
    //    if(bmAddIntervalSpans(fp, "1", starts+6, 20, values+6, 3)) goto error;
        //We can continue appending similarly formatted entries
    //    if(bmAppendIntervalSpans(fp, starts+9, values+9, 3)) goto error;

        //Add a new block of fixed-step entries
    //    if(bmAddIntervalSpanSteps(fp, "1", 900, 20, 30, values+12, 3)) goto error;
        //The start is then 760, since that's where the previous step ended
    //    if(bmAppendIntervalSpanSteps(fp, values+15, 3)) goto error;

        //Add a new chromosome
    //    chromsUse[0] = "2";
    //    chromsUse[1] = "2";
    //    chromsUse[2] = "2";
    //    if(bmAddIntervals(fp, chromsUse, starts, ends, values, 3)) goto error;

        //Closing the file causes the zoom levels to be created
        //if(DEBUG>0) 
        if(DEBUG>0) fprintf(stderr, "bm close1111 ----- \n");
        bmClose(fp);
        if(DEBUG>0) fprintf(stderr, "bm close22222 ---===--- \n");
        bmCleanup();
        /*
        * free memory from malloc
        */
    
        for(i =0; i < MAX_LINE_PRINT; i++){
            free(chroms[i]); free(chromsUse[i]); free(entryid[i]);
        }
        free(chroms); free(chromsUse); free(entryid); free(starts);
        free(ends); free(values); free(coverages); free(strands); free(contexts);
        free(chrom); free(old_chrom); free(chrLens); 
        free(strand); free(context); free(PerLine);free(chromlenf); 
        if(outfile) free(outfile); if(outbmfile) free(outbmfile);
        return 0;
    }

    /*
    * bigmeth file type
    * coverage: 1, strand: 2, context: 3, ID: 4
    * 0000000000001111, 0xf
    * 0000000011110000, 0xf0
    * 0000111100000000, 0xf00
    * 1111000000000000, 0xf000
    * int type; 0: bigmeth for loci [chrom start value coverage strand context]
    * int type; 1: bigmeth for loci [chrom start value coverage]
    * int type; 2: bigmeth for loci [chrom start value coverage context]
    * int type; 3: bigmeth for loci [chrom start value coverage strand]
    * 4: bigmeth for loci/region [chrom start end value coverage strand context]
    * 5: bigmeth for region with ID [chrom start end value coverage geneID]
    * 6: bigmeth for region without ID [chrom start end value coverage]
    * 7: bigmeth for region with ID and strand [chrom start end value coverage strand context geneID]
    */

    // read bm file, view
    if(strcmp(mode, "view")==0 || strcmp(mode, "addzm")==0 || strcmp(mode, "stats")==0){
        if(strcmp(mode, "addzm")==0){
            strcpy(outformat, "dm");
        }else if( strcmp(mode, "stats")==0 ){
            strcpy(outformat, "stats");
            memset(Fcover, 0, sizeof(int)*16);
        }
        if(DEBUG>0) fprintf(stderr, "dm view\n");
        //char region[] = "chr1:0-100,chr1:16766-16830";
        uint32_t type = BMtype(inbmfile, NULL);
        binaMethFile_t *ifp = NULL;
        ifp = bmOpen(inbmfile, NULL, "r");
        ifp->type = ifp->hdr->version;
        FILE *outfp_bm = NULL;
        FILE *outfp_stats = NULL;
        binaMethFile_t *ofp = NULL;
        if(strcmp(outformat, "txt") == 0 ){
            if(outfile) {
                outfp_bm = File_Open(outfile,"w");
            }else {
                outfp_bm = stdout;
            }
        }else if(strcmp(outformat, "stats") == 0){
            if(outfile) {
                outfp_bm = File_Open(outfile,"w");
                char *statsfile = malloc(200);
                strcpy(statsfile, outfile);
                strcat(statsfile, ".stats");
                outfp_stats = File_Open(statsfile,"w");
                free(statsfile);
            }else {
                outfp_bm = stdout;
                outfp_stats = stdout;
            }
        }else if(strcmp(outformat, "dm") == 0){
            fprintf(stderr, "output dm file\n");
            //bmfileinit(ofp, ifp, outfile, zoomlevel); // why not valid????
            if(bmInit(1<<17) != 0) {
                fprintf(stderr, "Received an error in bmInit\n");
                return 0;
            }
            ofp = bmOpen(outfile, NULL, "w");
            ofp->type = ifp->type;
            if(!ofp) {
                fprintf(stderr, "An error occurred while opening bm for writingn\n");
                return 0;
            }
            //Allow up to 10 zoom levels, though fewer will be used in practice
            if(bmCreateHdr(ofp, zoomlevel)) {
                fprintf(stderr, "== bmCreateHdr ==\n");
                return 0;
            }
            //Create the chromosome lists
            ofp->cl = bmCreateChromList_ifp(ifp); //2
            if(!ofp->cl) {
                fprintf(stderr, "== bmCreateChromList ==\n");
                return 0;
            }
            //Write the header
            if(bmWriteHdr(ofp)) {
                fprintf(stderr, "== bmWriteHdr ==\n");
                return 0;
            }
        }
        
        if(!region && !bedfile){
            fprintf(stderr, "[Mode] ------------- View all meth\n");
            main_view_all(ifp, outfp_bm, outformat, ofp, filterchrom);
        }else if(region){
            fprintf(stderr, "[Mode] ------------- View region %s\n", region);
            main_view(ifp, region, outfp_bm, outformat, ofp);
            free(region);
        }else if(bedfile){
            fprintf(stderr, "[Mode] ------------- View bedfile %s\n", bedfile);
            main_view_file(ifp, bedfile, outfp_bm, outformat, ofp);
            free(bedfile);
        }else{
            fprintf(stderr, "\nplease provide -r or --bed!!!\n");
        }

        if(strcmp(outformat, "stats") == 0) {
            int i = 0, j = 0;
            if(totalCG == 0){
                for(i=0; i<16; i++) totalCG += Fcover[i];
            }
            if(totalCG == 0) totalCG = 1;
            //fprintf(stderr, "%d\n", totalCG);
            for(i=0; i<16; i++){
                for(j=i+1; j<16; j++){
                    Fcover[i] += Fcover[j];
                }
            }

            fprintf(outfp_bm, "cover\t1");
            for(i=1; i<16; i++){
              fprintf(outfp_bm, "\t%d", i+1);
            }
            fprintf(outfp_bm, "\n");

            fprintf(outfp_bm, "count\t%d", Fcover[0]);
            for(i=1; i<16; i++){
              fprintf(outfp_bm, "\t%d", Fcover[i]);
            }
            fprintf(outfp_bm, "\n");

            fprintf(outfp_bm, "percent\t%.2f", (double)Fcover[0]/totalCG);
            for(i=1; i<16; i++){
              fprintf(outfp_bm, "\t%.2f", (double)Fcover[i]/totalCG);
            }
            fprintf(outfp_bm, "\n");

            // print stats
            //fprintf(outfp_stats, "m1\tm2\tm3\tm4\tm5\n%ld\t%ld\t%ld\t%ld\t%ld\n", mP1, mP2, mP3, mP4, mP5);
            fprintf(outfp_stats, "catary\t0");
            for(i=1; i<statsSize; i++){
              fprintf(outfp_stats, "\t%.2f", ((float)(i))/statsSize);
            }
            fprintf(outfp_stats, "\n");

            fprintf(outfp_stats, "percent\t%ld", mPs[0]);
            for(i=1; i< statsSize; i++){
                fprintf(outfp_stats, "\t%ld", mPs[i]);
            }
            fprintf(outfp_stats, "\n");
        }

        fprintf(stderr, "Done and free mem\n");
        free(inbmfile);
        bmClose(ifp);
        if(outfile) free(outfile);
        if(filterchrom) free(filterchrom);
        if(strcmp(outformat, "txt") == 0){
            fclose(outfp_bm);
        }else if(strcmp(outformat, "stats") == 0){
            fclose(outfp_bm);
            fclose(outfp_stats);
        }else if(strcmp(outformat, "dm") == 0){
            bmClose(ofp);
            bmCleanup();
        }
        return 0;
    }

    if(strcmp(mode, "viewheader")==0){
        if(DEBUG>0) fprintf(stderr, "dm viewheader\n");
        uint32_t type = BMtype(inbmfile, NULL);
        binaMethFile_t *ifp = NULL;
        ifp = bmOpen(inbmfile, NULL, "r");
        ifp->type = ifp->hdr->version;
        bmPrintHdr(ifp);
        //bmPrintIndexTree(ifp);
        free(inbmfile);
        bmClose(ifp);
        return 0;
    }

    //overlap
    if(strcmp(mode, "overlap")==0){
        fprintf(stderr, "[Mode] ------------- overlap %s %s\n", inbmfile, bmfile2);
        if(inbm_mul == 1){
            if(!region && !bedfile){
                bm_overlap_all_mul(inbmfiles, filter_strand);
            }else if(region){
                bm_overlap_region_mul(inbmfiles, region, filter_strand);
                free(region);
            }else if(bedfile){
                bm_overlap_file_mul(inbmfiles, bedfile);
                free(bedfile);
            }else{
                fprintf(stderr, "\nplease provide -r or --bed!!!\n");
            }
            free(inbmfiles); 
        }else{
            if(!region && !bedfile){
                bm_overlap_all(inbmfile, bmfile2, 1, 1, filter_strand);
            }else if(region){
                bm_overlap_region(inbmfile, bmfile2, region, filter_strand);
                free(region);
            }else if(bedfile){
                bm_overlap_file(inbmfile, bmfile2, bedfile);
                free(bedfile);
            }else{
                fprintf(stderr, "\nplease provide -r or --bed!!!\n");
            }
            free(inbmfile); 
            free(bmfile2); 
        }
        if(outfile) free(outfile);
        return 0;
    }

    // region
    if(strcmp(mode, "regionstats")==0){
        FILE *outfp_mean;
        if(outfile) {
            char* outfile_aver_temp = malloc(sizeof(char)*200);
            strcpy(outfile_aver_temp, outfile); strcat(outfile_aver_temp, ".regionm");
            outfp_mean = File_Open(outfile_aver_temp,"w");
            free(outfile_aver_temp);
        }else {
            outfp_mean = stdout;
        }

        FILE *outfp,*outfp_cg, *outfp_chg, *outfp_chh;
        if(printcoverage>0) {
            if(outfile){
                strcat(outfile , ".regionm.cover");
                char* outfile_temp = malloc(sizeof(char)*200);
                if(print2one != 0){ // dont write to one file
                    strcpy(outfile_temp, outfile);
                    outfp = File_Open(outfile_temp,"w");
                    outfp_cg = outfp;
                    outfp_chg = outfp;
                    outfp_chh = outfp;
                }else{
                    if(filter_context == 0 || filter_context == 4) {
                        strcpy(outfile_temp, outfile); strcat(outfile_temp, ".c");
                        outfp = File_Open(outfile_temp,"w");
                    }
                    if(filter_context == 0 || filter_context == 4 || filter_context == 1) {
                        strcpy(outfile_temp, outfile); strcat(outfile_temp, ".cg");
                        outfp_cg = File_Open(outfile_temp,"w");
                    }
                    if(filter_context == 0 || filter_context == 4 || filter_context == 2) {
                        strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chg");
                        outfp_chg = File_Open(outfile_temp,"w");
                    }
                    if(filter_context == 0 || filter_context == 4 || filter_context == 3) {
                        strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chh");
                        outfp_chh = File_Open(outfile_temp,"w");
                    }
                    free(outfile_temp);
                }
            }else{
                print2one = 2;
                outfp = stdout;
                outfp_cg = stdout;
                outfp_chg = stdout;
                outfp_chh = stdout;
            }
        }

        if(region){
            fprintf(stderr, "[Mode] ------------- regionstats %s %s %s\n", inbmfile, region, method);
            calregionstats(inbmfile, method, region, filter_strand, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(region);
        }else if(bedfile){
            fprintf(stderr, "[Mode] ------------- regionstats %s %s %s\n", inbmfile, bedfile, method);
            calregionstats_file(inbmfile, method, bedfile, 0, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(bedfile);
        }else if(gtffile){
            fprintf(stderr, "[Mode] ------------- regionstats %s %s %s\n", inbmfile, gtffile, method);
            calregionstats_file(inbmfile, method, gtffile, 1, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(gtffile);
        }else if(gfffile){
            fprintf(stderr, "[Mode] ------------- regionstats %s %s %s\n", inbmfile, gfffile, method);
            calregionstats_file(inbmfile, method, gfffile, 2, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(gfffile);
        }else{
            fprintf(stderr, "\nplease provide -r, --bed, --gtf or --gff!!!\n");
        }
        free(inbmfile);
        if(printcoverage>0) {
            if(filter_context == 0 || filter_context == 4) {
                fclose(outfp);
            }
            if(outfile && print2one == 0) { 
                if(filter_context == 0 || filter_context == 4 || filter_context == 1) fclose(outfp_cg);
                if(filter_context == 0 || filter_context == 4 || filter_context == 2) fclose(outfp_chg); 
                if(filter_context == 0 || filter_context == 4 || filter_context == 3) fclose(outfp_chh); 
            }
        }
        fclose(outfp_mean);
        if(outfile) free(outfile);
        return 0;
    }

    // region
    if(strcmp(mode, "bodystats")==0){
        FILE *outfp_mean;
        if(outfile) {
            char* outfile_aver_temp = malloc(sizeof(char)*200);
            strcpy(outfile_aver_temp, outfile); strcat(outfile_aver_temp, ".bodym");
            outfp_mean = File_Open(outfile_aver_temp,"w");
            free(outfile_aver_temp);
        }else {
            outfp_mean = stdout;
        }

        FILE *outfp, *outfp_cg, *outfp_chg, *outfp_chh;
        if(printcoverage>0) {
            if(outfile){
                strcat(outfile , ".bodym.cover");
                char* outfile_temp = malloc(sizeof(char)*200);
                if(print2one != 0){ // dont write to one file
                    strcpy(outfile_temp, outfile);
                    outfp = File_Open(outfile_temp,"w");
                    outfp_cg = outfp;
                    outfp_chg = outfp;
                    outfp_chh = outfp;
                }else{
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".c");
                    outfp = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".cg");
                    outfp_cg = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chg");
                    outfp_chg = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chh");
                    outfp_chh = File_Open(outfile_temp,"w");
                    free(outfile_temp);
                }
            }else{
                print2one = 2;
                outfp = stdout;
                outfp_cg = stdout;
                outfp_chg = stdout;
                outfp_chh = stdout;
            }
        }

        if(region){
            fprintf(stderr, "[Mode] ------------- bodystats %s %s %s\n", inbmfile, region, method);
            calbodystats(inbmfile, method, region, filter_strand, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(region);
        }else if(bedfile){
            fprintf(stderr, "[Mode] ------------- bodystats %s %s %s\n", inbmfile, bedfile, method);
            calbodystats_file(inbmfile, method, bedfile, 0, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(bedfile);
        }else if(gtffile){
            fprintf(stderr, "[Mode] ------------- bodystats %s %s %s\n", inbmfile, gtffile, method);
            calbodystats_file(inbmfile, method, gtffile, 1, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(gtffile);
        }else if(gfffile){
            fprintf(stderr, "[Mode] ------------- bodystats %s %s %s\n", inbmfile, gfffile, method);
            calbodystats_file(inbmfile, method, gfffile, 2, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
            free(gfffile);
        }else{
            fprintf(stderr, "\nplease provide -r, --bed, --gtf or --gff!!!\n");
        }

        free(inbmfile);
        if(printcoverage>0) {
            fclose(outfp);
            if(outfile && print2one == 0) { fclose(outfp_cg); fclose(outfp_chg); fclose(outfp_chh); }
        }
        fclose(outfp_mean);
        if(outfile) free(outfile);
        return 0;
    }

    // chromstats
    if(strcmp(mode, "chromstats")==0){
        fprintf(stderr, "[Mode] ------------- chromstats %s %s\n", inbmfile, method);
        FILE *outfp_mean;
        if(printcoverage == 0) {
            if(outfile) {
                char* outfile_aver_temp = malloc(sizeof(char)*200);
                strcpy(outfile_aver_temp, outfile); strcat(outfile_aver_temp, ".chrom");
                outfp_mean = File_Open(outfile_aver_temp,"w");
                free(outfile_aver_temp);
            }else {
                outfp_mean = stdout;
            }
        }

        FILE *outfp, *outfp_cg, *outfp_chg, *outfp_chh;
        if(printcoverage>0) {
            if(outfile){
                strcat(outfile , ".chrom.cover");
                char* outfile_temp = malloc(sizeof(char)*200);
                if(print2one != 0){ // dont write to one file
                    strcpy(outfile_temp, outfile);
                    outfp = File_Open(outfile_temp,"w");
                    outfp_cg = outfp;
                    outfp_chg = outfp;
                    outfp_chh = outfp;
                }else{
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".c");
                    outfp = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".cg");
                    outfp_cg = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chg");
                    outfp_chg = File_Open(outfile_temp,"w");
                    strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chh");
                    outfp_chh = File_Open(outfile_temp,"w");
                    free(outfile_temp);
                }
            }else{
                print2one = 2;
                outfp = stdout;
                outfp_cg = stdout;
                outfp_chg = stdout;
                outfp_chh = stdout;
            }
        }

        calchromstats(inbmfile, method, chromstep, stepoverlap, filter_strand, filter_context, outfp_mean, outfp, outfp_cg, outfp_chg, outfp_chh);
        free(inbmfile); 
        free(method);
        if(printcoverage>0) {
            fclose(outfp);
            if(outfile && print2one == 0) { fclose(outfp_cg); fclose(outfp_chg); fclose(outfp_chh); }
        }
        if(printcoverage == 0) fclose(outfp_mean); 
        if(outfile) free(outfile);
        return 0;
    }

    //gene profile
    if(strcmp(mode, "profile")==0){
        fprintf(stderr, "[Mode] ------------- profile %s\n", inbmfile);
        FILE *outfp, *outfp_cg, *outfp_chg, *outfp_chh;
        if(outfile){
            strcat(outfile ,profilemode_str); // add across, tss, center
            char* outfile_temp = malloc(sizeof(char)*200);
            if(print2one != 0){ // dont write to one file
                strcpy(outfile_temp, outfile);
                outfp = File_Open(outfile_temp,"w");
                outfp_cg = outfp;
                outfp_chg = outfp;
                outfp_chh = outfp;
            }else{
                strcpy(outfile_temp, outfile); strcat(outfile_temp, ".c");
                outfp = File_Open(outfile_temp,"w");
                strcpy(outfile_temp, outfile); strcat(outfile_temp, ".cg");
                outfp_cg = File_Open(outfile_temp,"w");
                strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chg");
                outfp_chg = File_Open(outfile_temp,"w");
                strcpy(outfile_temp, outfile); strcat(outfile_temp, ".chh");
                outfp_chh = File_Open(outfile_temp,"w");
                free(outfile_temp);
            }
        }else{
            print2one = 2;
            outfp = stdout;
            outfp_cg = stdout;
            outfp_chg = stdout;
            outfp_chh = stdout;
        }

        FILE *outfp_aver;
        if(outfile) {
            char* outfile_aver_temp = malloc(sizeof(char)*200);
            strcpy(outfile_aver_temp, outfile); strcat(outfile_aver_temp, ".aver");
            outfp_aver = File_Open(outfile_aver_temp,"w");
            free(outfile_aver_temp);
        }else {
            outfp_aver = stdout;
        }

        if(profilemode > 0) {
            profilestep = profilestep/2.0;
            profilemovestep = profilemovestep/2.0;
        }
        //bodyX
        bodyprofilestep = profilestep * bodyX;
        bodyprofilemovestep = profilemovestep * bodyX;

        if(bedfile){
            //calprofile(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, bedfile, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX, filter_strand);
            calprofile_gtf(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, bedfile, 3, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX, filter_strand);
            free(bedfile);
        }
        else if(gtffile){
            calprofile_gtf(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, gtffile, 1, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX, filter_strand);
            free(gtffile);
        }else if(gfffile){
            calprofile_gtf(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, gfffile, 2, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX, filter_strand);
            free(gfffile);
        }else{
            fprintf(stderr, "\nplease provide -r, --bed, --gtf or --gff!!!\n");
        }
        
        fclose(outfp); fclose(outfp_aver);
        if(outfile && print2one == 0) { fclose(outfp_cg); fclose(outfp_chg); fclose(outfp_chh); }
        if(outfile) free(outfile);
        return 0;
    }

    fprintf(stderr, "Please provide correct mode, unvalid mode: %s", mode);
    if(outfile) free(outfile);
    return 1;
error:
    fprintf(stderr, "Received an error in process!\n");
    bmClose(fp);
    bmCleanup();
    return -1;
}

void bmfileinit(binaMethFile_t *ofp, binaMethFile_t *ifp, char* outfile, int zoomlevel){
    if(bmInit(1<<17) != 0) {
        fprintf(stderr, "Received an error in bmInit\n");
        return;
    }
    ofp = bmOpen(outfile, NULL, "w");
    ofp->type = ifp->type;
    if(!ofp) {
        fprintf(stderr, "An error occurred while opening bm for writingn\n");
        return;
    }
    //Allow up to 10 zoom levels, though fewer will be used in practice
    if(bmCreateHdr(ofp, zoomlevel)) {
        fprintf(stderr, "== bmCreateHdr ==\n");
        return;
    }

    //Create the chromosome lists
    ofp->cl = ifp->cl; //2
    if(!ofp->cl) {
        fprintf(stderr, "== bmCreateChromList ==\n");
        return;
    }

    //Write the header
    if(bmWriteHdr(ofp)) {
        fprintf(stderr, "== bmWriteHdr ==\n");
        return;
    }
}

double *Sregionstats(binaMethFile_t *fp, char *chrom, int start, int end, int splitN, int binsize, uint32_t movestep, char *method, uint8_t strand, uint8_t context){
    double *stats = NULL;
    assert(splitN>0);
    //int i=0;
    if(strcmp(method, "mean")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, mean, strand, context);
    }else if(strcmp(method, "weighted")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, weighted, strand, context);
    }else if(strcmp(method, "dev")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, dev, strand, context);
    }else if(strcmp(method, "min")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, min, strand, context);
    }else if(strcmp(method, "max")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, max, strand, context);
    }else if(strcmp(method, "cover")==0){
        stats = bmStats(fp, chrom, start, end, splitN, binsize, movestep, cov, strand, context);
    }
    return stats;
}

//c cg chg chh
//idx0 c cg chg chh; idx1 c cg chg chh; ... etc
double *Sregionstats_array(binaMethFile_t *fp, char *chrom, int start, int end, int splitN, int binsize, uint32_t movestep, char *method, uint8_t strand){
    double *stats = NULL;
    assert(splitN>0);
    //int i=0;
    if(strcmp(method, "mean")==0){
        stats = bmStats_array(fp, chrom, start, end, splitN, binsize, movestep, mean, strand);
    }else if(strcmp(method, "weighted")==0){
        stats = bmStats_array(fp, chrom, start, end, splitN, binsize, movestep, weighted, strand);
    }
    return stats;
}

//double *output = malloc(sizeof(double)*nBins*Tsize);
void Sregionstats_array_count(binaMethFile_t *fp, char *chrom, int start, int end, int splitN, int binsize, uint32_t movestep, char *method, uint8_t strand, uint32_t *countC, uint32_t *countCT){
    assert(splitN>0);
    int i=0, Tsize = 4;
    for(i=0;i<splitN*Tsize;i++){
        countC[i]=0;
        countCT[i]=0;
    }
    if(strcmp(method, "weighted")==0){
        bmStats_array_count(fp, chrom, start, end, splitN, binsize, movestep, weighted, strand, countC, countCT);
    }else {
        fprintf(stderr, "Unexpected method for weighted methylation count!");
        exit(0);
    }
    return;
}

int calchromstats(char *inbmfile, char *method, int chromstep, int stepoverlap, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    int i = 0, start = 0, end = chromstep; //, j = 0
    char* region = malloc(sizeof(char)*1000);
    int splitN = 1, Tsize = 4;
    uint32_t *countC = malloc(sizeof(uint32_t)*splitN*Tsize);
    uint32_t *countCT = malloc(sizeof(uint32_t)*splitN*Tsize);
    for(i=0;i<fp->cl->nKeys;i++){
        char* chrom = (char*)fp->cl->chrom[i];
        int len = (int)fp->cl->len[i];
        start = 0; end = chromstep;
        //fprintf(stderr, "CCCC %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            //pstrand 0 +, 1 -, 2 ., 3 split
            if(pstrand>2){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, 0, 0, "", context, "+", outfileF);
                    calregion_print(fp, chrom, start, end, splitN, method, 1, 0, "", context, "-", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, 0, 0, "", context, "", "+", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, 1, 0, "", context, "", "-", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 0){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "+", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", "+", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 1){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "-", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", "-", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 2){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, ".", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", ".", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }
            start += stepoverlap;
            end += stepoverlap;
        }
    }
    free(region);

    bmClose(fp);
    bmCleanup();
    free(countC); free(countCT);
    return 0;
}

void profile_print(binaMethFile_t *fp, char *chrom, int start, int end, int splitN, int binsize, double profilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, FILE* outfileF){
    double *stats = Sregionstats(fp, chrom, start, end, splitN, binsize, (int)((end-start)*profilemovestep), "weighted", strand, context);
    int i;
    if(stats) {
        if(strand == 1) {// -
            fprintf(outfileF, "%s:%d-%d-%s\t%f", chrom, start, end, geneid, stats[splitN-1]);
            if(splitN>1){
                for(i=splitN-2;i>0;i--){
                    fprintf(outfileF, "\t%f", stats[i]);
                }
            }
            fprintf(outfileF, "\n");
        }else{
            fprintf(outfileF, "%s:%d-%d-%s\t%f", chrom, start, end, geneid, stats[0]);
            if(splitN>1){
                for(i=1;i<splitN;i++){
                    fprintf(outfileF, "\t%f", stats[i]);
                }
            }
            fprintf(outfileF, "\n");
        }
    }
}

void getcontext(uint8_t context, char* str_context){
    if(context == 0) strcpy(str_context, "C");
    else if(context == 1) strcpy(str_context, "CG");
    else if(context == 2) strcpy(str_context, "CHG");
    else if(context == 3) strcpy(str_context, "CHH");
    else strcpy(str_context, "Cx");
}

double *finalstats;
int *finalcounts;
FILE* outfileF;

typedef struct {
    binaMethFile_t *fp;
    char *chrom;
    int start;
    int end;
    int splitN;
    double profilemovestep;
    double bodyprofilemovestep;
    char *method;
    uint8_t strand;
    char* geneid;
    uint8_t context;
    int processtart;
    int processend;
    int bodysplitN;
    int matrixX;
    char* printbuffer_c;
    char* printbuffer_cg;
    char* printbuffer_chg;
    char* printbuffer_chh;
} ProfileARGS;

void profile_print_array(void *arg){
    int stored_buffer = 0;
    char* chrom = ((ProfileARGS*)arg)->chrom;
    int start = ((ProfileARGS*)arg)->start;
    int end = ((ProfileARGS*)arg)->end;
    int splitN = ((ProfileARGS*)arg)->splitN;
    double profilemovestep = ((ProfileARGS*)arg)->profilemovestep;
    char* method = ((ProfileARGS*)arg)->method;
    char* geneid = ((ProfileARGS*)arg)->geneid;
    uint8_t strand = ((ProfileARGS*)arg)->strand;
    uint8_t context = ((ProfileARGS*)arg)->context;
    int bodysplitN = ((ProfileARGS*)arg)->bodysplitN;
    int processtart = ((ProfileARGS*)arg)->processtart;
    int processend = ((ProfileARGS*)arg)->processend;
    double bodyprofilemovestep = ((ProfileARGS*)arg)->bodyprofilemovestep;
    int matrixX = ((ProfileARGS*)arg)->matrixX;
    char* printbuffer_c = ((ProfileARGS*)arg)->printbuffer_c;
    char* printbuffer_cg = ((ProfileARGS*)arg)->printbuffer_cg;
    char* printbuffer_chg = ((ProfileARGS*)arg)->printbuffer_chg;
    char* printbuffer_chh = ((ProfileARGS*)arg)->printbuffer_chh;

    //fprintf(stderr, "%s %d %d\n", chrom, start, end);
    int binsize = (int)((start-processtart)*profilemovestep)*2;
    double *stats_up = Sregionstats_array(((ProfileARGS*)arg)->fp, chrom, processtart, start, splitN, binsize, (int)((start-processtart)*profilemovestep), "weighted", strand);
    int bodybinsize = (int)((end-start)*bodyprofilemovestep)*2;
    double *stats_body = Sregionstats_array(((ProfileARGS*)arg)->fp, chrom, start, end, bodysplitN, bodybinsize, (int)((end-start)*bodyprofilemovestep), "weighted", strand);
    binsize = (int)((processend-end)*profilemovestep)*2;
    double *stats_down = Sregionstats_array(((ProfileARGS*)arg)->fp, chrom, end, processend, splitN, binsize, (int)((processend-end)*profilemovestep), "weighted", strand);
    int i,j,k,m; int Tsize = 4;
    int total_splitN_flank = splitN*Tsize;
    int total_splitN_body = bodysplitN*Tsize;
    //fprintf(stderr, "size %d %d\n", total_splitN_flank, total_splitN_body);
    int valid = 0;
    double *stats = malloc(sizeof(double)*(total_splitN_flank*2+total_splitN_body));
    if(stats_up && stats_body && stats_down){
        for(i=0;i<total_splitN_flank;i++){
            if(!isnan(stats_up[i])) {
                valid = 1;
                break;
            }
            if(!isnan(stats_down[i])) {
                valid = 1;
                break;
            }
        }
        for(i=0;i<total_splitN_body;i++){
            if(!isnan(stats_body[i])) {
                valid = 1;
                break;
            }
        }
    }
    if(stats_up && stats_body && stats_down){
        for(i=0;i<total_splitN_flank;i++){
            stats[i] = stats_up[i];
            stats[i+total_splitN_flank+total_splitN_body] = stats_down[i];
        }
        for(i=0;i<total_splitN_body;i++){
            stats[i+total_splitN_flank] = stats_body[i];
            //strtod("NaN", NULL);
        }
    }
    
    if(valid == 0) return;
    stored_buffer++;
    int total_splitN = total_splitN_flank*2+total_splitN_body;
    char* print_context = malloc(sizeof(char)*5);
    char* storetemp = malloc(sizeof(char)*100);
    float meth = 0; int mcount = 0;
    if(stats) {
        if(strand == 1) {// -
            for(j=total_splitN-Tsize, k=0; j<total_splitN; j++,k++){
                if(!(context>=4 || k==context)){
                    continue;
                }
                //getcontext(k, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                //store for average mr
                if(splitN>1){
                    for(i=j;i>=Tsize*(matrixX-1);){ //j-Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d, %d, %d %d %d\n", strand, i, total_splitN, j-i+m*Tsize+k, i-m*Tsize);
                            if(!isnan(stats[i-m*Tsize])) {
                                meth+=stats[i-m*Tsize];
                                mcount++;

                                finalstats[j-i+m*Tsize+k] += stats[i-m*Tsize];
                                finalcounts[j-i+m*Tsize+k]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
                    if(i>=0) {
                        // i >= 0
                        meth = 0, mcount = 0;
                        for(;i>=0;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[j+k-i] += stats[i];
                                finalcounts[j+k-i]++;
                            }
                            i=i-Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }
                }
                sprintf(storetemp, "\n");
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }else{ //+
            //0 C
            for(j=0; j<Tsize; j++){
                if(!(context>=4 || j==context)){
                    continue;
                }
                //getcontext(j, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                if(splitN>1){
                    for(i=j;i<total_splitN-Tsize*(matrixX-1);){//j+Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d %d %d\n", strand, i+m*Tsize, m);
                            if(!isnan(stats[i+m*Tsize])) {
                                meth+=stats[i+m*Tsize];
                                mcount++;

                                finalstats[i+m*Tsize] += stats[i+m*Tsize];
                                finalcounts[i+m*Tsize]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
                    // i < total_splitN
                    if(i<total_splitN) {
                        meth = 0, mcount = 0;
                        for(;i<total_splitN;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[i] += stats[i];
                                finalcounts[i]++;
                            }
                            i=i+Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }
                }
                sprintf(storetemp, "\n");
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
    return;
}

void profile_print_array_mp(void *arg){
    int stored_buffer = 0;
    char* chrom = ((ProfileARGS*)arg)->chrom;
    //int start = ((ProfileARGS*)arg)->start;
    //int end = ((ProfileARGS*)arg)->end;
    int splitN = ((ProfileARGS*)arg)->splitN;
    double profilemovestep = ((ProfileARGS*)arg)->profilemovestep;
    char* method = ((ProfileARGS*)arg)->method;
    char* geneid = ((ProfileARGS*)arg)->geneid;
    uint8_t strand = ((ProfileARGS*)arg)->strand;
    uint8_t context = ((ProfileARGS*)arg)->context;
    int bodysplitN = ((ProfileARGS*)arg)->bodysplitN;
    int processtart = ((ProfileARGS*)arg)->processtart;
    int processend = ((ProfileARGS*)arg)->processend;
    double bodyprofilemovestep = ((ProfileARGS*)arg)->bodyprofilemovestep;
    int matrixX = ((ProfileARGS*)arg)->matrixX;
    char* printbuffer_c = ((ProfileARGS*)arg)->printbuffer_c;
    char* printbuffer_cg = ((ProfileARGS*)arg)->printbuffer_cg;
    char* printbuffer_chg = ((ProfileARGS*)arg)->printbuffer_chg;
    char* printbuffer_chh = ((ProfileARGS*)arg)->printbuffer_chh;

    int binsize = (int)((processend-processtart)*profilemovestep)*2;
    double *stats = Sregionstats_array(((ProfileARGS*)arg)->fp, chrom, processtart, processend, splitN, binsize, (int)((processend-processtart)*profilemovestep), "weighted", strand);
    int i,j,k,m; int Tsize = 4;
    int total_splitN = splitN*Tsize;
    int valid = 0;
    if(stats){
        for(i=0;i<total_splitN;i++){
            if(!isnan(stats[i])) {
                valid = 1;
                break;
            }
        }
    }

    if(valid == 0) return;
    stored_buffer++;
    char* print_context = malloc(sizeof(char)*5);
    char* storetemp = malloc(sizeof(char)*100);
    float meth = 0; int mcount = 0;
    if(stats) {
        if(strand == 1) {// -
            for(j=total_splitN-Tsize, k=0; j<total_splitN; j++,k++){
                if(!(context>=4 || k==context)){
                    continue;
                }
                //getcontext(k, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, processtart, processend, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, processtart, processend);
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                //store for average mr
                if(splitN>1){
                    for(i=j;i>=Tsize*(matrixX-1);){ //j-Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d, %d, %d %d %d\n", strand, i, total_splitN, j-i+m*Tsize+k, i-m*Tsize);
                            if(!isnan(stats[i-m*Tsize])) {
                                meth+=stats[i-m*Tsize];
                                mcount++;

                                finalstats[j-i+m*Tsize+k] += stats[i-m*Tsize];
                                finalcounts[j-i+m*Tsize+k]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
                    if(i>=0) {
                        // i >= 0
                        meth = 0, mcount = 0;
                        for(;i>=0;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[j+k-i] += stats[i];
                                finalcounts[j+k-i]++;
                            }
                            i=i-Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }

                }
                sprintf(storetemp, "\n");
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }else{ //+
            //0 C
            for(j=0; j<Tsize; j++){
                if(!(context>=4 || j==context)){
                    continue;
                }
                //getcontext(j, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, processtart, processend, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, processtart, processend);
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                if(splitN>1){
                    for(i=j;i<total_splitN-Tsize*(matrixX-1);){//j+Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d %d %d\n", strand, i+m*Tsize, m);
                            if(!isnan(stats[i+m*Tsize])) {
                                meth+=stats[i+m*Tsize];
                                mcount++;

                                finalstats[i+m*Tsize] += stats[i+m*Tsize];
                                finalcounts[i+m*Tsize]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
                    if(i<total_splitN) {
                        // i < total_splitN
                        meth = 0, mcount = 0;
                        for(;i<total_splitN;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[i] += stats[i];
                                finalcounts[i]++;
                            }
                            i=i+Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }

                }
                sprintf(storetemp, "\n");
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
}

uint32_t stored_buffer = 0;
void profile_print_array_buffer(double *finalstats, int *finalcounts, binaMethFile_t *fp, char *chrom, int processtart, int start, int end, int processend, int splitN, int bodysplitN, double profilemovestep, double bodyprofilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, char* printbuffer_c, char* printbuffer_cg, char* printbuffer_chg, char* printbuffer_chh, int matrixX){
    int binsize = (int)((start-processtart)*profilemovestep)*2;
    double *stats_up = Sregionstats_array(fp, chrom, processtart, start, splitN, binsize, (int)((start-processtart)*profilemovestep), "weighted", strand);
    binsize = (int)((end-start)*bodyprofilemovestep)*2;
    double *stats_body = Sregionstats_array(fp, chrom, start, end, bodysplitN, binsize, (int)((end-start)*bodyprofilemovestep), "weighted", strand);
    binsize = (int)((processend-end)*profilemovestep)*2;
    double *stats_down = Sregionstats_array(fp, chrom, end, processend, splitN, binsize, (int)((processend-end)*profilemovestep), "weighted", strand);
    int i,j,k,m; int Tsize = 4;
    int total_splitN_flank = splitN*Tsize;
    int total_splitN_body = bodysplitN*Tsize;
    //fprintf(stderr, "size %d %d\n", total_splitN_flank, total_splitN_body);
    int valid = 0;
    double *stats = malloc(sizeof(double)*(total_splitN_flank*2+total_splitN_body));
    if(stats_up && stats_body && stats_down){
        for(i=0;i<total_splitN_flank;i++){
            if(!isnan(stats_up[i])) {
                valid = 1;
                break;
            }
            if(!isnan(stats_down[i])) {
                valid = 1;
                break;
            }
        }
        for(i=0;i<total_splitN_body;i++){
            if(!isnan(stats_body[i])) {
                valid = 1;
                break;
            }
        }
    }
    if(stats_up && stats_body && stats_down){
        for(i=0;i<total_splitN_flank;i++){
            stats[i] = stats_up[i];
            stats[i+total_splitN_flank+total_splitN_body] = stats_down[i];
        }
        for(i=0;i<total_splitN_body;i++){
            stats[i+total_splitN_flank] = stats_body[i];
            //strtod("NaN", NULL);
        }
    }
    
    if(valid == 0) return;
    stored_buffer++;
    int total_splitN = total_splitN_flank*2+total_splitN_body;
    char* print_context = malloc(sizeof(char)*5);
    char* storetemp = malloc(sizeof(char)*100);
    float meth = 0; int mcount = 0;
    if(stats) {
        if(strand == 1) {// -
            for(j=total_splitN-Tsize, k=0; j<total_splitN; j++,k++){
                if(!(context>=4 || k==context)){
                    continue;
                }
                //getcontext(k, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                //store for average mr
                if(splitN>1){
                    for(i=j;i>=Tsize*(matrixX-1);){ //j-Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d, %d, %d %d %d\n", strand, i, total_splitN, j-i+m*Tsize+k, i-m*Tsize);
                            if(!isnan(stats[i-m*Tsize])) {
                                meth+=stats[i-m*Tsize];
                                mcount++;

                                finalstats[j-i+m*Tsize+k] += stats[i-m*Tsize];
                                finalcounts[j-i+m*Tsize+k]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
                    if(i>=0) {
                        // i >= 0
                        meth = 0, mcount = 0;
                        for(;i>=0;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[j+k-i] += stats[i];
                                finalcounts[j+k-i]++;
                            }
                            i=i-Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }
                }
                sprintf(storetemp, "\n");
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }else{ //+
            //0 C
            for(j=0; j<Tsize; j++){
                if(!(context>=4 || j==context)){
                    continue;
                }
                //getcontext(j, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                if(splitN>1){
                    for(i=j;i<total_splitN-Tsize*(matrixX-1);){//j+Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d %d %d\n", strand, i+m*Tsize, m);
                            if(!isnan(stats[i+m*Tsize])) {
                                meth+=stats[i+m*Tsize];
                                mcount++;

                                finalstats[i+m*Tsize] += stats[i+m*Tsize];
                                finalcounts[i+m*Tsize]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
                    if(i<total_splitN) {
                        // i < total_splitN
                        meth = 0, mcount = 0;
                        for(;i<total_splitN;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[i] += stats[i];
                                finalcounts[i]++;
                            }
                            i=i+Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }
                }
                sprintf(storetemp, "\n");
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
}

void profile_print_array_buffer1(double *finalstats, int *finalcounts, binaMethFile_t *fp, char *chrom, int start, int end, int splitN, double profilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, char* printbuffer_c, char* printbuffer_cg, char* printbuffer_chg, char* printbuffer_chh, int matrixX){
    int binsize = (int)((end-start)*profilemovestep)*2;
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, binsize, (int)((end-start)*profilemovestep), "weighted", strand);
    int i,j,k,m; int Tsize = 4;
    int total_splitN = splitN*Tsize;
    int valid = 0;
    if(stats){
        for(i=0;i<total_splitN;i++){
            if(!isnan(stats[i])) {
                valid = 1;
                break;
            }
        }
    }

    if(valid == 0) return;
    stored_buffer++;
    char* print_context = malloc(sizeof(char)*5);
    char* storetemp = malloc(sizeof(char)*100);
    float meth = 0; int mcount = 0;
    if(stats) {
        if(strand == 1) {// -
            for(j=total_splitN-Tsize, k=0; j<total_splitN; j++,k++){
                if(!(context>=4 || k==context)){
                    continue;
                }
                //getcontext(k, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                //store for average mr
                if(splitN>1){
                    for(i=j;i>=Tsize*(matrixX-1);){ //j-Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d, %d, %d %d %d\n", strand, i, total_splitN, j-i+m*Tsize+k, i-m*Tsize);
                            if(!isnan(stats[i-m*Tsize])) {
                                meth+=stats[i-m*Tsize];
                                mcount++;

                                finalstats[j-i+m*Tsize+k] += stats[i-m*Tsize];
                                finalcounts[j-i+m*Tsize+k]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
                    if(i>=0) {
                        // i >= 0
                        meth = 0, mcount = 0;
                        for(;i>=0;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[j+k-i] += stats[i];
                                finalcounts[j+k-i]++;
                            }
                            i=i-Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }

                }
                sprintf(storetemp, "\n");
                if(k==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(k==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(k==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(k==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }else{ //+
            //0 C
            for(j=0; j<Tsize; j++){
                if(!(context>=4 || j==context)){
                    continue;
                }
                //getcontext(j, print_context);
                if(geneid[0]) sprintf(storetemp, "%s:%d-%d-%s", chrom, start, end, geneid);
                else sprintf(storetemp, "%s:%d-%d", chrom, start, end);
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                if(splitN>1){
                    for(i=j;i<total_splitN-Tsize*(matrixX-1);){//j+Tsize
                        meth = 0, mcount = 0;
                        for(m=0;m<matrixX;m++){
                            //fprintf(stderr, "%d %d %d\n", strand, i+m*Tsize, m);
                            if(!isnan(stats[i+m*Tsize])) {
                                meth+=stats[i+m*Tsize];
                                mcount++;

                                finalstats[i+m*Tsize] += stats[i+m*Tsize];
                                finalcounts[i+m*Tsize]++;
                            }
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
                    if(i<total_splitN) {
                        // i < total_splitN
                        meth = 0, mcount = 0;
                        for(;i<total_splitN;){
                            if(!isnan(stats[i])) {
                                meth+=stats[i];
                                mcount++;

                                finalstats[i] += stats[i];
                                finalcounts[i]++;
                            }
                            i=i+Tsize;
                        }
                        if(mcount>0) sprintf(storetemp, "\t%.4lf", meth/mcount);
                        else strcpy(storetemp, "\tnan");
                        if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                        else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                        else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                        else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
                    }

                }
                sprintf(storetemp, "\n");
                if(j==0) printbuffer_c = fastStrcat(printbuffer_c, storetemp);
                else if(j==1) printbuffer_cg = fastStrcat(printbuffer_cg, storetemp);
                else if(j==2) printbuffer_chg = fastStrcat(printbuffer_chg, storetemp);
                else if(j==3) printbuffer_chh = fastStrcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
}

int calprofile_gtf(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *gtffile, int format, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX, uint8_t filt_strand){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fgtffile=File_Open(gtffile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(200*sizeof(char));
    char *strand = malloc(2); 
    int pstrand = 2; //.
    char *geneid = malloc(200*sizeof(char));
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilemovestep)-1; //*3
    int bodysplitN = ceil(1.0/bodyprofilemovestep)-1;
    assert(splitN>0 && bodysplitN>0);

    int Total_splitN = 0;
    if(profilemode == 0) {
        Total_splitN = splitN*2 + bodysplitN;
    }else{
        Total_splitN = splitN; //*2;
    }
    //aver
    int Tsize = 4;
    finalstats = calloc(Total_splitN*Tsize, sizeof(double));
    finalcounts = calloc(Total_splitN*Tsize, sizeof(int));
    int start=0, end=0, middle = 0, processtart = 0, processend = 0;
    //char* storebuffer_c = malloc(sizeof(char)*MAX_BUFF_PRINT);
    //char* storebuffer_cg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    //char* storebuffer_chg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    //char* storebuffer_chh = malloc(sizeof(char)*MAX_BUFF_PRINT);
    uint32_t Nprocess = 0;
    //NTHREAD = 10
    int prothreads = 0, ithread = 0;
    ProfileARGS *profileArgs = malloc(sizeof(ProfileARGS)*NTHREAD);
    //binaMethFile_t **fps = malloc(sizeof(binaMethFile_t*)*NTHREAD);
    for(ithread=0; ithread<NTHREAD; ithread++){
        if(profilemode>0){
            profileArgs[ithread].splitN = splitN;//*2;
        }else{
            profileArgs[ithread].splitN = splitN;
        }
        profileArgs[ithread].bodysplitN = bodysplitN;
        profileArgs[ithread].bodyprofilemovestep = bodyprofilemovestep;
        profileArgs[ithread].profilemovestep = profilemovestep;
        profileArgs[ithread].matrixX = matrixX;
        profileArgs[ithread].method = "weighted";
        profileArgs[ithread].fp = NULL;
        profileArgs[ithread].fp = bmOpen(inbmfile, NULL, "r");
        profileArgs[ithread].fp->type = profileArgs[ithread].fp->hdr->version;
        profileArgs[ithread].printbuffer_c = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_cg = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_chg = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_chh = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].geneid = malloc(200*sizeof(char));
        profileArgs[ithread].chrom = malloc(200*sizeof(char));
    }
    pthread_t* Threads=(pthread_t*) malloc(sizeof(pthread_t)*NTHREAD);

    while(fgets(PerLine,2000,Fgtffile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 3){ // bed
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
            strcpy(geneid, "");
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
        }
        //fprintf(stderr, "%s %d %d %s %s\n", chrom, start, end, strand, geneid);
        if(end-start < bodysplitN) continue;
        Nprocess++;
        if(Nprocess%5000==0) fprintf(stderr, "Processed %d regions\n", Nprocess);

        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
            fprintf(stderr, "Waring: undetected strand infor\n");
        }
        if(filt_strand!=2 && filt_strand!=pstrand && pstrand!=2) continue;

        strcpy(profileArgs[prothreads].chrom, chrom);
        profileArgs[prothreads].start = start;
        profileArgs[prothreads].end = end;
        profileArgs[prothreads].strand = pstrand;
        strcpy(profileArgs[prothreads].geneid, geneid);
        profileArgs[prothreads].context = context;
        if(profilemode == 0){ // gene flanks mode
            if(start > upstream) processtart = start-upstream;
            else processtart = 0;
            processend = end+downstream;
            profileArgs[prothreads].processtart = processtart;
            profileArgs[prothreads].processend = processend;
            
            //fprintf(stderr, "TT %s %d %d %d\n", chrom, start, end, prothreads);
            //profile_print_array(&profileArgs);
            int r=pthread_create(&Threads[prothreads],NULL,profile_print_array,(void*) &profileArgs[prothreads]);
            if(r) {fprintf(stderr, "Launch_Threads():Cannot create thread..\n");exit(-1);}
            
            prothreads++;
            if(prothreads>=NTHREAD){
                for (ithread=0;ithread<NTHREAD;ithread++)
                {
                    pthread_join(Threads[ithread],NULL);
                }
                prothreads = 0;
            }
            //profile_print_array_buffer(finalstats, finalcounts, fp, chrom, processtart, start, end, processend, splitN, bodysplitN, profilemovestep, bodyprofilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 1) { // gene tss mode
            if(strand[0] == '+' || strand[0] == '.') {
                if(start > upstream) processtart = start- upstream;
                else processtart = 0;
                processend = start+downstream;
            }else{
                if(end > upstream) processtart = end - upstream;
                else processtart = 0;
                processend = end+downstream;
            }
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 2) { // gene tts mode
            if(strand[0] == '+' || strand[0] == '.') {
                if(end > upstream) processtart = end - upstream;
                else processtart = 0;
                processend = end+downstream;
            }else{
                if(start > upstream) processtart = start- upstream;
                else processtart = 0;
                processend = start+downstream;
            }
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 3){ //center mode
            middle = (int)((start+end)/2);
            if(middle>upstream) processtart = middle - upstream;
            else processtart = 0;
            processend = middle + downstream;
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }
        if(profilemode>0){
            profileArgs[prothreads].processtart = processtart;
            profileArgs[prothreads].processend = processend;

            //fprintf(stderr, "TT %s %d %d %d\n", chrom, start, end, prothreads);
            int r=pthread_create(&Threads[prothreads],NULL,profile_print_array_mp,(void*) &profileArgs[prothreads]);
            if(r) {fprintf(stderr, "Launch_Threads():Cannot create thread..\n");exit(-1);}
            
            prothreads++;
            if(prothreads>=NTHREAD){
                for (ithread=0;ithread<NTHREAD;ithread++)
                {
                    pthread_join(Threads[ithread],NULL);
                }
                prothreads = 0;
            }
        }
        // print result
        for (ithread=0;ithread<NTHREAD;ithread++)
        {
            fprintf(outfileF_c,"%s",profileArgs[ithread].printbuffer_c);
            profileArgs[ithread].printbuffer_c[0] = '\0';
            fprintf(outfileF_cg,"%s",profileArgs[ithread].printbuffer_cg);
            profileArgs[ithread].printbuffer_cg[0] = '\0';
            fprintf(outfileF_chg,"%s",profileArgs[ithread].printbuffer_chg);
            profileArgs[ithread].printbuffer_chg[0] = '\0';
            fprintf(outfileF_chh,"%s",profileArgs[ithread].printbuffer_chh);
            profileArgs[ithread].printbuffer_chh[0] = '\0';
            stored_buffer = 0;
        }
    }

    fprintf(stderr, "print aver matrix\n");
    //print aver
    char* print_context = malloc(sizeof(char)*5);
    int j;
    for(j=0; j<Tsize; j++){
        if(!(context>=4 || j==context)){
            continue;
        }
        getcontext(j, print_context);
        if(finalcounts[j]>0) fprintf(outfileF_aver, "%s\t%f", print_context, finalstats[j]/finalcounts[j]);
        else fprintf(outfileF_aver, "%s\t%f", print_context, strtod("NaN", NULL));
        if(Total_splitN>1){
            for(i=j+Tsize;i<Total_splitN*Tsize;){
                if(finalcounts[i]>0) fprintf(outfileF_aver, "\t%f", finalstats[i]/finalcounts[i]);
                else fprintf(outfileF_aver, "\t%f", strtod("NaN", NULL));
                i+=Tsize;
            }
        }
        fprintf(outfileF_aver, "\n");
    }
    fprintf(stderr, "free and close file\n");
    free(print_context);
    free(finalstats); free(finalcounts);//???
    fclose(Fgtffile);
    free(chrom); free(PerLine); free(strand); free(geneid); 
    bmClose(fp);
    for(ithread=0; ithread<NTHREAD; ithread++){
        //bmClose(profileArgs[ithread].fp);
        free(profileArgs[ithread].printbuffer_c);
        free(profileArgs[ithread].printbuffer_cg);
        free(profileArgs[ithread].printbuffer_chg);
        free(profileArgs[ithread].printbuffer_chh);
        free(profileArgs[ithread].geneid);
        free(profileArgs[ithread].chrom);
    }
    free(profileArgs);
    bmCleanup();
    fprintf(stderr, "done free!\n");
}

int calprofile(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *bedfile, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX, uint8_t filt_strand){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(200*sizeof(char));
    char *strand = malloc(2); int pstrand = 2; //.
    char *geneid = malloc(200*sizeof(char));
    strcpy(geneid, "");
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilemovestep)-1; //*3
    int bodysplitN = ceil(1.0/bodyprofilemovestep)-1;
    assert(splitN>0 && bodysplitN>0);
    int Total_splitN = 0;
    if(profilemode == 0) {
        Total_splitN = splitN*2 + bodysplitN;
    }else{
        Total_splitN = splitN;//*2;
    }
    //aver
    int Tsize = 4;
    finalstats = calloc(Total_splitN*Tsize, sizeof(double));
    finalcounts = calloc(Total_splitN*Tsize, sizeof(int));
    int start=0, end=0, middle = 0, processtart = 0, processend = 0;
    uint16_t Nprocess = 0;

    //int NTHREAD = 10
    int prothreads = 0, ithread = 0;
    ProfileARGS *profileArgs = malloc(sizeof(ProfileARGS)*NTHREAD);
    //binaMethFile_t **fps = malloc(sizeof(binaMethFile_t*)*NTHREAD);
    for(ithread=0; ithread<NTHREAD; ithread++){
        if(profilemode>0){
            profileArgs[ithread].splitN = splitN;//*2;
        }else{
            profileArgs[ithread].splitN = splitN;
        }
        profileArgs[ithread].bodysplitN = bodysplitN;
        profileArgs[ithread].bodyprofilemovestep = bodyprofilemovestep;
        profileArgs[ithread].profilemovestep = profilemovestep;
        profileArgs[ithread].matrixX = matrixX;
        profileArgs[ithread].method = "weighted";
        profileArgs[ithread].fp = NULL;
        profileArgs[ithread].fp = bmOpen(inbmfile, NULL, "r");
        profileArgs[ithread].fp->type = profileArgs[ithread].fp->hdr->version;
        profileArgs[ithread].printbuffer_c = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_cg = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_chg = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].printbuffer_chh = malloc(sizeof(char)*MAX_BUFF_PRINT);
        profileArgs[ithread].chrom = malloc(200*sizeof(char));
        profileArgs[ithread].geneid = malloc(200*sizeof(char));
    }
    pthread_t* Threads=(pthread_t*) malloc(sizeof(pthread_t)*NTHREAD);

    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
            strand[0] = '.';
            strand[1] = '\0';
        }
        if(end-start < splitN) continue;
        Nprocess++;
        if(Nprocess%5000==0) fprintf(stderr, "Processed %d regions\n", Nprocess);

        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        if(filt_strand!=2 && filt_strand!=pstrand && pstrand!=2) continue;

        strcpy(profileArgs[prothreads].chrom, chrom);
        profileArgs[prothreads].start = start;
        profileArgs[prothreads].end = end;
        profileArgs[prothreads].strand = pstrand;
        strcpy(profileArgs[prothreads].geneid, geneid);
        profileArgs[prothreads].context = context;

        if(profilemode == 0){ // gene flanks mode
            if(start > upstream) processtart = start-upstream;
            else processtart = 0;
            processend = end+downstream;
            
            profileArgs[prothreads].processtart = processtart;
            profileArgs[prothreads].processend = processend;
            
            //fprintf(stderr, "TT %s %d %d %d\n", chrom, start, end, prothreads);
            //profile_print_array(&profileArgs);
            int r=pthread_create(&Threads[prothreads],NULL,profile_print_array,(void*) &profileArgs[prothreads]);
            if(r) {fprintf(stderr, "Launch_Threads():Cannot create thread..\n");exit(-1);}
            
            prothreads++;
            if(prothreads>=NTHREAD){
                for (ithread=0;ithread<NTHREAD;ithread++)
                {
                    pthread_join(Threads[ithread],NULL);
                }
                prothreads = 0;
            }
            //profile_print_array_buffer(finalstats, finalcounts, fp, chrom, processtart, start, end, processend, splitN, bodysplitN, profilemovestep, bodyprofilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 1) { // gene tss mode
            if(strand[0] == '+' || strand[0] == '.') {
                if(start > upstream) processtart = start- upstream;
                else processtart = 0;
                processend = start+downstream;
            }else{
                if(end > upstream) processtart = end - upstream;
                else processtart = 0;
                processend = end+downstream;
            }
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 2) { // gene tts mode
            if(strand[0] == '+' || strand[0] == '.') {
                if(end > upstream) processtart = end - upstream;
                else processtart = 0;
                processend = end+downstream;
            }else{
                if(start > upstream) processtart = start- upstream;
                else processtart = 0;
                processend = start+downstream;
            }
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }else if(profilemode == 3){ //center mode
            middle = (int)((start+end)/2);
            if(middle>upstream) processtart = middle - upstream;
            else processtart = 0;
            processend = middle + downstream;
            //profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
        }
        if(profilemode>0){
            profileArgs[prothreads].processtart = processtart;
            profileArgs[prothreads].processend = processend;

            //fprintf(stderr, "TT %s %d %d %d\n", chrom, start, end, prothreads);
            int r=pthread_create(&Threads[prothreads],NULL,profile_print_array_mp,(void*) &profileArgs[prothreads]);
            if(r) {fprintf(stderr, "Launch_Threads():Cannot create thread..\n");exit(-1);}
            
            prothreads++;
            if(prothreads>=NTHREAD){
                for (ithread=0;ithread<NTHREAD;ithread++)
                {
                    pthread_join(Threads[ithread],NULL);
                }
                prothreads = 0;
            }
        }
        // print result
        for (ithread=0;ithread<NTHREAD;ithread++)
        {
            fprintf(outfileF_c,"%s",profileArgs[ithread].printbuffer_c);
            profileArgs[ithread].printbuffer_c[0] = '\0';
            fprintf(outfileF_cg,"%s",profileArgs[ithread].printbuffer_cg);
            profileArgs[ithread].printbuffer_cg[0] = '\0';
            fprintf(outfileF_chg,"%s",profileArgs[ithread].printbuffer_chg);
            profileArgs[ithread].printbuffer_chg[0] = '\0';
            fprintf(outfileF_chh,"%s",profileArgs[ithread].printbuffer_chh);
            profileArgs[ithread].printbuffer_chh[0] = '\0';
            stored_buffer = 0;
        }
    }

    fprintf(stderr, "print aver matrix\n");
    //print aver
    char* print_context = malloc(sizeof(char)*5);
    int j;
    for(j=0; j<Tsize; j++){
        if(!(context>=4 || j==context)){
            continue;
        }
        getcontext(j, print_context);
        if(finalcounts[j]>0) fprintf(outfileF_aver, "%s\t%f", print_context, finalstats[j]/finalcounts[j]);
        else fprintf(outfileF_aver, "%s\t%f", print_context, strtod("NaN", NULL));
        if(Total_splitN>1){
            for(i=j+Tsize;i<Total_splitN*Tsize;){
                if(finalcounts[i]>0) fprintf(outfileF_aver, "\t%f", finalstats[i]/finalcounts[i]);
                else fprintf(outfileF_aver, "\t%f", strtod("NaN", NULL));
                i+=Tsize;
            }
        }
        fprintf(outfileF_aver, "\n");
    }

    fprintf(stderr, "free and close file\n");
    fclose(Fbedfile);
    free(print_context);
    free(finalstats); free(finalcounts);

    bmClose(fp);
    free(chrom); free(PerLine); free(strand); free(geneid);
    for(ithread=0; ithread<NTHREAD; ithread++){
        //bmClose(profileArgs[ithread].fp);
        //why core dump this four line... //printbuffer_cxx
        free(profileArgs[ithread].printbuffer_c);
        free(profileArgs[ithread].printbuffer_cg);
        free(profileArgs[ithread].printbuffer_chg);
        free(profileArgs[ithread].printbuffer_chh);
        free(profileArgs[ithread].chrom);
        //free(profileArgs[ithread].geneid);
    }
    free(profileArgs);
    bmCleanup();
}

void delete_char(char str[],char target){
	int i,j;
	for(i=j=0;str[i]!='\0';i++){
		if(str[i]!=target){
			str[j++]=str[i];
		}
	}
	str[j]='\0';
}

void delete_char2(char* str, char target,char target2){
	int i,j;
	for(i=j=0;str[i]!='\0';i++){
		if(str[i]!=target && str[i]!=target2 ){
			str[j++]=str[i];
		}
	}
	str[j]='\0';
}

void calbody_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF){
    int binsize = end-start;
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, binsize, end-start, method, pstrand);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    getcontext(context, print_context);
    if(stats) {
        if(format == 1 || format == 2){
            if(context >= 4){
                if(!isnan(stats[0])) {
                    fprintf(outfileF, "C\t%s\t%f\t%s\t%s\n", bodycase, stats[0], strand, geneid);
                }
                if(!isnan(stats[1])) {
                    fprintf(outfileF, "CG\t%s\t%f\t%s\t%s\n", bodycase, stats[1], strand, geneid);
                }
                if(!isnan(stats[2])) {
                    fprintf(outfileF, "CHG\t%s\t%f\t%s\t%s\n", bodycase, stats[2], strand, geneid);
                }
                if(!isnan(stats[3])) {
                    fprintf(outfileF, "CHH\t%s\t%f\t%s\t%s\n", bodycase, stats[3], strand, geneid);
                }
            }else if(!isnan(stats[context])) {
                fprintf(outfileF, "%s\t%s\t%f\t%s\t%s\n", print_context, bodycase, stats[context], strand, geneid);
            }
        }else{
            if(context >= 4) {
                if(!isnan(stats[0])) {
                    fprintf(outfileF, "C\t%s\t%f\t%s:%d-%d\t%s\n", bodycase, stats[0], chrom, start, end, strand);
                }
                if(!isnan(stats[1])) {
                    fprintf(outfileF, "CG\t%s\t%f\t%s:%d-%d\t%s\n", bodycase, stats[1], chrom, start, end, strand);
                }
                if(!isnan(stats[2])) {
                    fprintf(outfileF, "CHG\t%s\t%f\t%s:%d-%d\t%s\n", bodycase, stats[2], chrom, start, end, strand);
                }
                if(!isnan(stats[3])) {
                    fprintf(outfileF, "CHH\t%s\t%f\t%s:%d-%d\t%s\n", bodycase, stats[3], chrom, start, end, strand);
                }
            }else if(!isnan(stats[context])) {
                fprintf(outfileF, "%s\t%s\t%f\t%s:%d-%d\t%s\n", print_context, bodycase, stats[context], chrom, start, end, strand);
            }
        }
        /*
        if(splitN>1){
            for(i=1;i<splitN;i++){
                printf("\t%f", stats[i]);
            }
        }*/
    }
    free(print_context);
}

void calregion_weighted_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, uint32_t *countC, uint32_t *countCT){
    int binsize = end-start;
    Sregionstats_array_count(fp, chrom, start, end, splitN, binsize, end-start, method, pstrand, countC, countCT);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    char* print_geneid = malloc(sizeof(char)*200);
    if(format == 1 || format == 2) {
        strcpy(print_geneid, "\t");
        strcat(print_geneid, geneid);
    }else{
        sprintf(print_geneid, "\t%s:%d-%d", chrom, start, end);
    }
    getcontext(context, print_context);
    if(countCT && countC) {
        if(format == 1 || format == 2 || format == 0){ //bed gtf etc
            if(context >= 4){
                if(countCT[0]>0) { // C
                    fprintf(outfileF_c, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[0], countCT[0], strand[0], print_geneid);
                }
                if(countCT[1]>0) { // CG
                    fprintf(outfileF_cg, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[1], countCT[1], strand[0], print_geneid);
                }
                if(countCT[2]>0) {
                    fprintf(outfileF_chg, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[2], countCT[2], strand[0], print_geneid);
                }
                if(countCT[3]>0) {
                    fprintf(outfileF_chh, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[3], countCT[3], strand[0], print_geneid);
                }
            }else if(countCT[context]>0) { // print_context
                if(context == 0) fprintf(outfileF_c, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[context], countCT[context], strand[0], print_geneid);
                else if(context == 1) fprintf(outfileF_cg, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[context], countCT[context], strand[0], print_geneid);
                else if(context == 2) fprintf(outfileF_chg, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[context], countCT[context], strand[0], print_geneid);
                else if(context == 3) fprintf(outfileF_chh, "%s\t%d\t%d\t%d\t%d\t%c%s\n", chrom, start, end, countC[context], countCT[context], strand[0], print_geneid);
            }
        }

    }
    free(print_context);
    free(print_geneid);
}

void calregion_print(binaMethFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* strand, FILE* outfileF){
    int binsize = end-start;
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, binsize, end-start, method, pstrand);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    getcontext(context, print_context);
    if(stats) {
        if(format == 1 || format == 2){ // gtf gff
            if(context >= 4){
                if(!isnan(stats[0])) { //C
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tC\t%c\t%s\n", chrom, start, end, stats[0], strand[0], geneid);
                }
                if(!isnan(stats[1])) { //CG
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCG\t%c\t%s\n", chrom, start, end, stats[1], strand[0], geneid);
                }
                if(!isnan(stats[2])) { //CHG
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCHG\t%c\t%s\n", chrom, start, end, stats[2], strand[0], geneid);
                }
                if(!isnan(stats[3])) { //CHH
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCHH\t%c\t%s\n", chrom, start, end, stats[3], strand[0], geneid);
                }
            }else if(!isnan(stats[context])) { //print_context
                fprintf(outfileF, "%s\t%d\t%d\t%f\t%s\t%c\t%s\n", chrom, start, end, stats[context], print_context, strand[0], geneid);
            }else if(alwaysprint == 1){
                fprintf(outfileF, "%s\t%d\t%d\tNA\t%s\t%c\t%s\n", chrom, start, end, print_context, strand[0], geneid);
            }
        }else{ //bed or region
            if(context >= 4){
                if(!isnan(stats[0])) { //C
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tC\t%c\n", chrom, start, end, stats[0], strand[0]);
                }
                if(!isnan(stats[1])) { //CG
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCG\t%c\n", chrom, start, end, stats[1], strand[0]);
                }
                if(!isnan(stats[2])) { //CHG
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCHG\t%c\n", chrom, start, end, stats[2], strand[0]);
                }
                if(!isnan(stats[3])) { //CHH
                    fprintf(outfileF, "%s\t%d\t%d\t%f\tCHH\t%c\n", chrom, start, end, stats[3], strand[0]);
                }
            }else if(!isnan(stats[context])) { //print_context
                fprintf(outfileF, "%s\t%d\t%d\t%f\t%s\t%c\n", chrom, start, end, stats[context], print_context, strand[0]);
            }else if(alwaysprint == 1){
                fprintf(outfileF, "%s\t%d\t%d\tNA\t%s\t%c\n", chrom, start, end, print_context, strand[0]);
            }
        }
    }
    free(print_context);
}

int calregionstats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, Tsize = 4; //, i =0
    char *geneid = malloc(200*sizeof(char));
    uint32_t *countC = malloc(sizeof(uint32_t)*splitN*Tsize);
    uint32_t *countCT = malloc(sizeof(uint32_t)*splitN*Tsize);
    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 0){
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
            if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
                strand[0] = '.';
                strand[1] = '\0';
            }
        }else if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
        }
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        //fprintf(stderr, "%s\n", PerLine);
        //fprintf(stderr, "%s%d%d%s", chrom, start, end, strand);
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
    
        if(printcoverage == 0){
            calregion_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, strand, outfileF);
        }else{
            calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, "", strand, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
        }
    }

    bmClose(fp);
    bmCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
    free(countC); free(countCT);
    return 0;
}


int calregionstats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0; //, j = 0
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    int splitN = 1, Tsize = 4;
    char *strandstr = malloc(100*sizeof(char)); //int strand = 2;
    uint8_t my_pstrand = 2;
    uint16_t *countC = malloc(sizeof(uint16_t)*splitN*Tsize);
    uint16_t *countCT = malloc(sizeof(uint16_t)*splitN*Tsize);
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        strandstr = strtok(NULL,",:-");
        if(strandstr) {
            if(strandstr[0] == '+'){
                my_pstrand = 0;
            }else if(strandstr[0] == '-'){
                my_pstrand = 1;
            }else{
                my_pstrand = 2;
            }
        
            if(printcoverage == 0){
                calregion_print(fp, chrom, start, end, splitN, method, my_pstrand, 0, "", context, strandstr, outfileF);
            }else{
                calregion_weighted_print(fp, chrom, start, end, splitN, method, my_pstrand, 0, "", context, "", strandstr, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
            }
        }else{
            if(pstrand>2){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, 0, 0, "", context, "+", outfileF);
                    calregion_print(fp, chrom, start, end, splitN, method, 1, 0, "", context, "-", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, 0, 0, "", context, "", "+", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, 1, 0, "", context, "", "-", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 0){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "+", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", "+", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 1){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "-", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", "-", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }else if(pstrand == 2){
                if(printcoverage == 0){
                    calregion_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, ".", outfileF);
                }else{
                    calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", ".", outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
                }
            }
        }
        
    }

    bmClose(fp);
    bmCleanup();
    //free(chrom);
    free(countC); free(countCT);
    return 0;
}

int calbodystats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, Tsize = 4; //, i =0
    char *geneid = malloc(200*sizeof(char));
    int upstream = 0, downstream = 0;
    uint32_t *countC = malloc(sizeof(uint32_t)*splitN*Tsize);
    uint32_t *countCT = malloc(sizeof(uint32_t)*splitN*Tsize);
    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 0){ //bed
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
            if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
                strand[0] = '.';
                strand[1] = '\0';
            }
        }else if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){ //gff
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
        }
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        //fprintf(stderr, "%s\n", PerLine);
        //fprintf(stderr, "%s%d%d%s", chrom, start, end, strand);
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        upstream = start-regionextend>=0?start-regionextend:0;
        downstream = end + regionextend;
        calbody_print(fp, chrom, upstream, start, splitN, method, pstrand, format, geneid, context, "up", strand, outfileF);
        calbody_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, "body", strand, outfileF);
        calbody_print(fp, chrom, end, downstream, splitN, method, pstrand, format, geneid, context, "down", strand, outfileF);
        if(printcoverage > 0){
            calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, "body", strand, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
        }
    }

    bmClose(fp);
    bmCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
    free(countC); free(countCT);
    return 0;
}


int calbodystats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open bm file
    uint32_t type = BMtype(inbmfile, NULL);
    binaMethFile_t *fp = NULL;
    fp = bmOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0;//, j = 0
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    int splitN = 1, Tsize = 4;
    char *strandstr = malloc(100*sizeof(char)); //int strand = 2;
    int upstream = 0, downstream = 0;
    uint32_t *countC = malloc(sizeof(uint32_t)*splitN*Tsize);
    uint32_t *countCT = malloc(sizeof(uint32_t)*splitN*Tsize);
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        strandstr = strtok(NULL,",:-");
        if(strandstr) {
            if(strandstr[0] == '+'){
                pstrand = 0;
            }else if(strandstr[0] == '-'){
                pstrand = 1;
            }else{
                pstrand = 2;
            }
        }
        upstream = start-regionextend>=0?start - regionextend:0;
        downstream = end + regionextend;
        calbody_print(fp, chrom, upstream, start, splitN, method, pstrand, 0, "", context, "up", strandstr, outfileF);
        calbody_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "body", strandstr, outfileF);
        calbody_print(fp, chrom, end, downstream, splitN, method, pstrand, 0, "", context, "down", strandstr, outfileF);
    
        if(printcoverage > 0){
            calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "body", strandstr, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
        }
    }

    bmClose(fp);
    bmCleanup();
    //free(chrom);
    free(countC); free(countCT);
    return 0;
}

void write_dm(binaMethFile_t *ifp, char* region, FILE* outfileF, char *outformat, binaMethFile_t *ofp, char **chromsUse, uint32_t *starts, uint32_t *ends, float *values, uint16_t *coverages, uint8_t *strands, uint8_t *contexts, char **entryid, int newchr){
    int printL = 0;
    printL = main_view_bm(ifp, region, outfileF, outformat, ofp, chromsUse, starts, ends, values, coverages, strands, contexts,
         entryid);
    if(strcmp(outformat, "dm") == 0) {
        //if(start == 0 && printL>0){
        //这里需要加个判断，要写入的区域是否已经存在相同染色体，如果存在用bmAppendIntervals
        if(newchr == 0 && printL>0){
            int response = bmAddIntervals(ofp, chromsUse, starts, ends, values, coverages, strands, contexts,
            entryid, printL);
            if(response) {
                fprintf(stderr, "bmAddIntervals 0\n");
                return -1;
            }
            printL = 0;
        }else if(printL>0){
            if(bmAppendIntervals(ofp, starts, ends, values, coverages, strands, contexts, entryid, printL)) {
                fprintf(stderr, "bmAppendIntervals 2\n");
                return -1;
            }
            printL = 0;
        }
    }
}

int main_view_all(binaMethFile_t *ifp, FILE* outfileF, char *outformat, binaMethFile_t *ofp, char* filterchrom){

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    int i=0, printL = 0;
    char* region = malloc(sizeof(char)*1000);

    char **chromsUse = malloc(sizeof(char*)*MAX_LINE_PRINT);
    char **entryid = malloc(sizeof(char*)*MAX_LINE_PRINT);
    uint32_t *chrLens = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    uint32_t *starts = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    uint32_t *ends = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
    float *values = malloc(sizeof(float) * MAX_LINE_PRINT);
    uint16_t *coverages = malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
    uint8_t *strands = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
    uint8_t *contexts = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);

    for(i=0;i<ifp->cl->nKeys;i++){
        char* chrom = (char*)ifp->cl->chrom[i];
        int len = (int)ifp->cl->len[i];
        start = 0, end = SEGlen-1;
        if(filterchrom){
            if(strcmp(chrom, filterchrom)!=0){
                continue;
            }
        }
        fprintf(stderr, "process %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            sprintf(region, "%s:%d-%d\n", chrom, start, end);
            //fprintf(stderr, "ccx %s\n", region);
            write_dm(ifp, region, outfileF, outformat, ofp, chromsUse, starts, ends, values, coverages, strands, contexts, entryid, start);

            //printL = main_view_bm(ifp, region, outfileF, outformat, ofp, chromsUse, starts, ends, values, coverages, strands, contexts, 
            //    entryid);
            //if(strcmp(outformat, "dm") == 0) {
            //    if(start == 0 && printL>0){
            //        int response = bmAddIntervals(ofp, chromsUse, starts, ends, values, coverages, strands, contexts, 
            //        entryid, printL);
            //        if(response) {
            //            fprintf(stderr, "bmAddIntervals 0\n");
            //            return -1;
            //        }
            //        printL = 0;
            //    }else if(printL>0){
            //        if(bmAppendIntervals(ofp, starts, ends, values, coverages, strands, contexts, entryid, printL)) {
            //            fprintf(stderr, "bmAppendIntervals 2\n");
            //            return -1;
            //        }
            //        printL = 0;
            //    }
            //}
            start += SEGlen;
            end += SEGlen;
        }
    }
//    fprintf(stderr, "\nWWW0W\n");

    free(region);
    for(i =0; i < MAX_LINE_PRINT; i++){
        free(chromsUse[i]); free(entryid[i]);
    }
    free(chromsUse); free(entryid); free(starts);
    free(ends); free(values); free(coverages); free(strands); free(contexts);
}

int main_view_bm(binaMethFile_t *ifp, char *region, FILE* outfileF, char *outformat, binaMethFile_t *ofp, \
    char** chromsUse, uint32_t* starts, uint32_t* ends, float* values, uint16_t* coverages, uint8_t* strands, \
    uint8_t* contexts, char** entryid){
    // read. test/example_output.bm
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", region, ifp->type);
    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *pszBuf = malloc(sizeof(char)*40000000);
    char *tempstore = pszBuf;
    char *tempchar = malloc(30);
    int Nprint = 0; int cover = 0; float methlevel = 0;
    //char *strand = malloc(100*sizeof(char)); int strand;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        //strand = strtok(NULL,",:-");
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        bmOverlappingIntervals_t *o;
        o = bmGetOverlappingIntervals(ifp, chrom, start, end+1);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        //fprintf(stderr, "--- version --- %d %d\n", ifp->hdr->version, o->l);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                if(ifp->hdr->version & BM_COVER){
                    if(!(o->coverage[j]>=mincover && o->coverage[j]<=maxcover)){
                        continue;
                    }
                }
                if(filter_strand != 2 && ifp->hdr->version & BM_STRAND) {
                    if(o->strand[j] != filter_strand){
                        continue;
                    }
                }
                if(filter_context!=0 && filter_context<4 && ifp->hdr->version & BM_CONTEXT){
                    if(o->context[j] != filter_context){
                        continue;
                    }
                }

                if(strcmp(outformat, "dm") == 0){
                    chromsUse[Nprint] = strdup(chrom);
                    starts[Nprint] = o->start[j];
                    if(ifp->hdr->version & BM_END){
                        ends[Nprint] = o->end[j];
                    }else{
                        ends[Nprint] = starts[Nprint] + 1;
                    }
                    values[Nprint] = o->value[j];
                    if(ifp->hdr->version & BM_COVER) {
                        coverages[Nprint] = o->coverage[j];
                    }
                    if(ifp->hdr->version & BM_STRAND){
                        strands[Nprint] = o->strand[j];
                    }
                    if(ifp->hdr->version & BM_CONTEXT) {
                        contexts[Nprint] = o->context[j];
                    }
                    if(ifp->hdr->version & BM_ID) {
                        strcpy(entryid[Nprint], o->entryid[j]);
                    }
                    Nprint++;
                }else if(strcmp(outformat, "stats") == 0) {
                    methlevel = o->value[j];
                    //if(methlevel>=0.8) mP5++;
                    //else if(methlevel >= 0.6) mP4++;
                    //else if(methlevel >= 0.4) mP3++;
                    //else if(methlevel >= 0.2) mP2++;
                    //else mP1++;
                    mPs[(int)(methlevel*statsSize - 0.01/statsSize)]++;

                    if(ifp->hdr->version & BM_COVER) {
                        cover = o->coverage[j]-1;
                        if(cover>=0 && cover < 15){
                            Fcover[cover]++;
                        }else Fcover[15]++;
                    }else{
                        fprintf(stderr, "Undetected coverage information in DM file, please check and reprint!\n");
                        exit(0);
                    }
                }else if(strcmp(outformat, "txt") == 0){
                    //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                    //o->strand[i], o->context[i]);   //%"PRIu32"
                    sprintf(tempchar, "%s\t%ld", chrom, o->start[j]);
                    tempstore = fastStrcat(tempstore, tempchar);
                    if(ifp->hdr->version & BM_END) { //ifp->type
                        sprintf(tempchar, "\t%"PRIu32"", o->end[j]);
                        tempstore = fastStrcat(tempstore, tempchar);
                    }
                    sprintf(tempchar, "\t%f", o->value[j]);
                    tempstore = fastStrcat(tempstore, tempchar);
                    if(ifp->hdr->version & BM_COVER) {
                        sprintf(tempchar, "\t%"PRIu16"", o->coverage[j]);
                        tempstore = fastStrcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_STRAND){
                        sprintf(tempchar, "\t%s", strand_str[o->strand[j]]);
                        tempstore = fastStrcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_CONTEXT) {
                        sprintf(tempchar, "\t%s", context_str[o->context[j]]);
                        tempstore = fastStrcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_ID) {
                        sprintf(tempchar, "\t%s", o->entryid[j]);
                        tempstore = fastStrcat(tempstore, tempchar);
                    }
                    sprintf(tempchar, "\n");
                    tempstore = fastStrcat(tempstore, tempchar);
                    Nprint++;
                    if(Nprint>10000) {
                        fprintf(outfileF,"%s",pszBuf);
                        Nprint = 0;
                        pszBuf[0] = '\0';
                        tempstore = pszBuf;
                    }
                }
            }
        }
        bmDestroyOverlappingIntervals(o);
    }

    if(Nprint>0) {
        if(strcmp(outformat, "txt") == 0) {
            fprintf(outfileF,"%s",pszBuf);
        } else if(strcmp(outformat, "dm") == 0) {
            ;//bm out, just skip
        }
    }
    //free(chrom);
    free(tempchar);
    free(pszBuf);
    return Nprint;

error:
    fprintf(stderr, "No results found!\n");
    return 0;
}

int main_view_file(binaMethFile_t *ifp, char *bedfile, FILE* outfileF, char *outformat, binaMethFile_t *ofp){
    // read. test/example_output.bm
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;
    bmOverlappingIntervals_t *o;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", bedfile, ifp->type);

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    unsigned int j = 0;
    char *tempstore = malloc(sizeof(char)*10000000);
    char *tempchar = malloc(20);
    int Nprint = 0;
    char* region = malloc(sizeof(char)*1000);
    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
             strand[0] = '.';
             strand[1] = '\0';
        }
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen chrom %s %d %d",  chrom, start, end);
        sprintf(region, "%s:%d-%d", chrom, start, end);
        //fprintf(stderr, "--- version --- %d\n", ifp->hdr->version);
        main_view(ifp, region, outfileF, outformat, ofp);
    }

    fclose(Fbedfile);
    free(chrom);
    free(strand);
    free(PerLine); free(region);

    bmDestroyOverlappingIntervals(o);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    return 1;
}

int main_view(binaMethFile_t *ifp, char *region, FILE* outfileF, char *outformat, binaMethFile_t *ofp){
    // read. test/example_output.bm
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", region, ifp->type);
    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *pszBuf = malloc(sizeof(char)*40000000);
    char *tempstore = pszBuf;
    char *tempchar = malloc(20);
    int Nprint = 0; int cover = 0; float methlevel = 0;
    if(strcmp(outformat, "dm") == 0) {
        char **chromsUse = malloc(sizeof(char*)*MAX_LINE_PRINT);
        char **entryid = malloc(sizeof(char*)*MAX_LINE_PRINT);
        uint32_t *chrLens = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
        uint32_t *starts = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
        uint32_t *ends = malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
        float *values = malloc(sizeof(float) * MAX_LINE_PRINT);
        uint16_t *coverages = malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
        uint8_t *strands = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
        uint8_t *contexts = malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
        for(i=0;i<slen; i++){
            //chrom = strtok(regions[i], ",:-");
            //start = atoi(strtok(NULL,",:-"));
            //end = atoi(strtok(NULL,",:-"));
            if(end-start>1000000) {
                fprintf(stderr, "view region bigger than 1Mb, exit;");
                exit(0);
            }
            write_dm(ifp, regions[i], outfileF, outformat, ofp, chromsUse, starts, ends, values, coverages, strands, contexts, entryid, 0);
        }
        //free mem
        for(i =0; i < MAX_LINE_PRINT; i++){
            free(chromsUse[i]); free(entryid[i]);
        }
        free(chromsUse); free(entryid); free(starts);
        free(ends); free(values); free(coverages); free(strands); free(contexts);
    }else{
        //char *strand = malloc(100*sizeof(char)); int strand;
        for(i=0;i<slen; i++){
            chrom = strtok(regions[i], ",:-");
            start = atoi(strtok(NULL,",:-"));
            end = atoi(strtok(NULL,",:-")); // + 1;
            if(end-start>1000000) {
                fprintf(stderr, "view region bigger than 1Mb, exit;");
                exit(0);
            }
            //strand = strtok(NULL,",:-");
            //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
            if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);

            bmOverlappingIntervals_t *o;
            o = bmGetOverlappingIntervals(ifp, chrom, start, end+1);
            if(!o) goto error;
            if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
            //fprintf(stderr, "--- version --- %d\n", ifp->hdr->version);
            if(o->l) {
                for(j=0; j<o->l; j++) {
                    if(ifp->hdr->version & BM_COVER){
                        if(!(o->coverage[j]>=mincover && o->coverage[j]<=maxcover)){
                            continue;
                        }
                    }
                    if(filter_strand != 2 && ifp->hdr->version & BM_STRAND) {
                        if(o->strand[j] != filter_strand){
                            continue;
                        }
                    }
                    if(filter_context!=0 && filter_context<4 && ifp->hdr->version & BM_CONTEXT){
                        if(o->context[j] != filter_context){
                            continue;
                        }
                    }

                    if(strcmp(outformat, "stats") == 0) {
                        methlevel = o->value[j];
                        //if(methlevel>=0.8) mP5++;
                        //else if(methlevel >= 0.6) mP4++;
                        //else if(methlevel >= 0.4) mP3++;
                        //else if(methlevel >= 0.2) mP2++;
                        //else mP1++;
                        mPs[(int)(methlevel*statsSize - 0.01/statsSize)]++;

                        if(ifp->hdr->version & BM_COVER) {
                            cover = o->coverage[j]-1;
                            if(cover>=0 && cover < 15){
                                Fcover[cover]++;
                            }else Fcover[15]++;
                        }else{
                            fprintf(stderr, "Undetected coverage information in DM file, please check and reprint!\n");
                            exit(0);
                        }
                    }else if(strcmp(outformat, "txt") == 0) {
                        //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                        //o->strand[i], o->context[i]);   %"PRIu32"
                        sprintf(tempchar, "%s\t%ld", chrom, o->start[j]);
                        tempstore = fastStrcat(tempstore, tempchar);
                        if(ifp->hdr->version & BM_END) { //ifp->type
                            sprintf(tempchar, "\t%"PRIu32"", o->end[j]);
                            tempstore = fastStrcat(tempstore, tempchar);
                        }
                        sprintf(tempchar, "\t%f", o->value[j]);
                        tempstore = fastStrcat(tempstore, tempchar);
                        if(ifp->hdr->version & BM_COVER) {
                            sprintf(tempchar, "\t%"PRIu16"", o->coverage[j]);
                            tempstore = fastStrcat(tempstore, tempchar);
                        }
                        if(ifp->hdr->version & BM_STRAND){
                            sprintf(tempchar, "\t%s", strand_str[o->strand[j]]);
                            tempstore = fastStrcat(tempstore, tempchar);
                        }
                        if(ifp->hdr->version & BM_CONTEXT) {
                            sprintf(tempchar, "\t%s", context_str[o->context[j]]);
                            tempstore = fastStrcat(tempstore, tempchar);
                        }
                        if(ifp->hdr->version & BM_ID) {
                            sprintf(tempchar, "\t%s", o->entryid[j]);
                            tempstore = fastStrcat(tempstore, tempchar);
                        }
                        sprintf(tempchar, "\n");
                        tempstore = fastStrcat(tempstore, tempchar);
                        Nprint++;
                        if(Nprint>10000) {
                            fprintf(outfileF,"%s",pszBuf);
                            Nprint = 0;
                            pszBuf[0] = '\0';
                            tempstore = pszBuf;
                        }
                    }
                }
            }
            bmDestroyOverlappingIntervals(o);
        }
    }

    if(Nprint>0) {
        if(strcmp(outformat, "txt") == 0) fprintf(outfileF,"%s",pszBuf);
        else if(strcmp(outformat, "dm") == 0) {
            //bm out
        }
        Nprint = 0;
    }
    //free(chrom);
    free(tempchar);
    free(pszBuf);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    return 1;
}

int main_view_bedfile(char *inbmF, char *bedfile, int type, FILE* outfileF, char *outformat, binaMethFile_t *ofp){
    // read. test/example_output.bm
    binaMethFile_t *ifp = NULL;
    ifp = bmOpen(inbmF, NULL, "r");
    ifp->type = type;
    bmOverlappingIntervals_t *o;

    FILE *BedF = File_Open(bedfile, "r");

    if(DEBUG>1) fprintf(stderr, "xxx1211 %s\n", bedfile);
    int i =0, j = 0;
    
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char region[200]; int cover = 0; float methlevel = 0;
    while (fgets(region,200,BedF)!=0){
        chrom = strtok(region, ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-"));
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d chrom %s %d %d", i, chrom, start, end);
        o = bmGetOverlappingIntervals(ifp, chrom, start, end+1);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                if(strcmp(outformat, "stats") == 0) {
                    methlevel = o->value[j];
                    //if(methlevel>=0.8) mP5++;
                    //else if(methlevel >= 0.6) mP4++;
                    //else if(methlevel >= 0.4) mP3++;
                    //else if(methlevel >= 0.2) mP2++;
                    //else mP1++;
                    mPs[(int)(methlevel*statsSize - 0.01/statsSize)]++;

                    if(ifp->hdr->version & BM_COVER) {
                        cover = o->coverage[j]-1;
                        if(cover>=0 && cover < 15){
                            Fcover[cover]++;
                        }else Fcover[15]++;
                    }else{
                        fprintf(stderr, "Undetected coverage information in DM file, please check and reprint!\n");
                        exit(0);
                    }
                }else if(strcmp(outformat, "txt") == 0) {
                    //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                    //o->strand[i], o->context[i]);
                    fprintf(stderr, "%s\t%ld\t%ld\t%f\t%ld\t%s\t%s\t%s\n", chrom, o->start[j], o->end[j], o->value[j], o->coverage[j],
                        strand_str[o->strand[j]], context_str[o->context[j]], o->entryid[j]);
                }
            }
        }
    }

    //free(chrom);
    fclose(BedF);
    bmDestroyOverlappingIntervals(o);
    bmClose(ifp);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    bmClose(ifp);
    bmCleanup();
    return 1;
}

int bm_overlap_all_mul(char *inbmFs, uint8_t pstrand){
    char *substr= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        substr = strtok(NULL,",");
    }

    binaMethFile_t **ifps = malloc(sizeof(binaMethFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        binaMethFile_t *ifp1 = NULL;
        ifp1 = bmOpen(infiles[i], NULL, "r");
        ifp1->type = ifp1->hdr->version;
        ifps[i] = ifp1;
    }

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    for(i=0;i<ifps[0]->cl->nKeys;i++){
        char* chrom = (char*)ifps[0]->cl->chrom[i];
        int len = (int)ifps[0]->cl->len[i];
        start = 0, end = SEGlen-1;
        //fprintf(stderr, "CCCC %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            bm_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
            start += SEGlen;
            end += SEGlen;
        }
    }

    
    for(i=0;i<sizeifp;i++){
        bmClose(ifps[i]);
    }
    return 0;
}

int bm_overlap_file_mul(char *inbmFs, char *bedfile){
    /*
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    binaMethFile_t *ifp1 = NULL;
    ifp1 = bmOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    binaMethFile_t *ifp2 = NULL;
    ifp2 = bmOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;
    */


    char *substr= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        substr = strtok(NULL,",");
    }

    binaMethFile_t **ifps = malloc(sizeof(binaMethFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        binaMethFile_t *ifp1 = NULL;
        ifp1 = bmOpen(infiles[i], NULL, "r");
        ifp1->type = ifp1->hdr->version;
        ifps[i] = ifp1;
    }
    
    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int SEGlen = 1000000, nK = 0, j = 0;
    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
            strand[0] = '.';
            strand[1] = '\0';
        }
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        if(end-start>SEGlen){
            nK = (int)((end-start)/SEGlen-0.5);
            for(j=0;j<nK;j++){
                bm_overlap_mul(ifps, sizeifp, chrom, start+j*SEGlen, start+(j+1)*SEGlen-1, pstrand);
            }
            bm_overlap_mul(ifps, sizeifp, chrom, start+nK*SEGlen, end, pstrand);
        }else
            bm_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
    }

    for(i=0;i<sizeifp;i++){
        bmClose(ifps[i]);
    }
    fclose(Fbedfile);
    free(PerLine); free(chrom); free(strand);
}

int bm_overlap_file(char *inbmF1, char *inbmF2, char *bedfile){
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    binaMethFile_t *ifp1 = NULL;
    ifp1 = bmOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    binaMethFile_t *ifp2 = NULL;
    ifp2 = bmOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    
    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(2000);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    while(fgets(PerLine,2000,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(!(strand[0]=='+' || strand[0]=='-' || strand[0]=='.')) {
            strand[0] = '.';
            strand[1] = '\0';
        }
        if(end<start){
            fprintf(stderr, "Warning, chromosome start bigger than end, %s %d %d", chrom, start, end);
            continue;
        }
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        bm_overlap(ifp1, ifp2, chrom, start, end, pstrand);
    }

    bmClose(ifp1);
    bmClose(ifp2);
    fclose(Fbedfile);
    free(PerLine); free(chrom); free(strand);
}

int bm_overlap_region_mul(char *inbmFs, char *region, uint8_t pstrand){
    //open file1
    char *subinbm= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (subinbm != NULL) {
        strcpy(infiles[sizeifp++], subinbm);
        subinbm = strtok(NULL,",");
    }

    binaMethFile_t **ifps = malloc(sizeof(binaMethFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        binaMethFile_t *ifp1 = NULL;
        ifp1 = bmOpen(infiles[i], NULL, "r");
        ifp1->type = ifp1->hdr->version;
        ifps[i] = ifp1;
    }

    //region
    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    int SEGlen = 1000000;
    char *chrom = malloc(100*sizeof(char));
    int start=0, end=0, nK = 0, j = 0;
    char *strandstr = malloc(100*sizeof(char));
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        strandstr = strtok(NULL,",:-");
        if(strandstr) {
            if(strandstr[0] == '+'){
                pstrand = 0;
            }else if(strandstr[0] == '-'){
                pstrand = 1;
            }else{
                pstrand = 2;
            }
        }
        if(end-start>SEGlen){
            nK = (int)((end-start)/SEGlen-0.5);
            for(j=0;j<nK;j++){
                bm_overlap_mul(ifps, sizeifp, chrom, start+j*SEGlen, start+(j+1)*SEGlen-1, pstrand);
            }
            bm_overlap_mul(ifps, sizeifp, chrom, start+nK*SEGlen, end, pstrand);
        }else
            bm_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
    }

    for(i=0;i<sizeifp;i++){
        bmClose(ifps[i]);
    }
    return 0;
}

int bm_overlap_region(char *inbmF1, char *inbmF2, char *region, uint8_t pstrand){
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    binaMethFile_t *ifp1 = NULL;
    ifp1 = bmOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    binaMethFile_t *ifp2 = NULL;
    ifp2 = bmOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    char *substr= strtok(region, ",");
    char regions[1000][200] = {""};
    int slen = 0, i =0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char));
    int start=0, end=0;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-")); // + 1;
        bm_overlap(ifp1, ifp2, chrom, start, end, pstrand);
    }

    bmClose(ifp1);
    bmClose(ifp2);
    return 0;
}

int bm_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand){
    if(n1<1 || n2<1){
        return -1;
    }
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    binaMethFile_t *ifp1 = NULL;
    ifp1 = bmOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    binaMethFile_t *ifp2 = NULL;
    ifp2 = bmOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    int i=0;
    for(i=0;i<ifp1->cl->nKeys;i++){
        char* chrom = (char*)ifp1->cl->chrom[i];
        int len = (int)ifp1->cl->len[i];
        start = 0, end = SEGlen-1;
        //fprintf(stderr, "CCCC %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            bm_overlap(ifp1, ifp2, chrom, start, end, pstrand);
            start += SEGlen;
            end += SEGlen;
        }
    }

    bmClose(ifp1);
    bmClose(ifp2);
    return 0;
}

int bm_overlap(binaMethFile_t *ifp1, binaMethFile_t *ifp2, char *chrom, int start, int end, uint8_t strand){

    int slen = 1, i =0, j = 0, k = 0, lociK = 0;
    int* countM = malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    for(i=0;i<slen; i++){
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        bmOverlappingIntervals_t *o1;
        bmOverlappingIntervals_t *o2;
        o1 = bmGetOverlappingIntervals(ifp1, chrom, start, end+1);
        o2 = bmGetOverlappingIntervals(ifp2, chrom, start, end+1);
        if(!o1 || !o2) goto error;
        //fprintf(stderr, "--- version --- %d\n", ifp1->hdr->version);
        if(o1->l && o2->l) {
            for(j=0; j<o1->l; j++) {
                if(strand!=2){
                    if(ifp1->hdr->version & BM_STRAND){
                        if(strand != o1->strand[j]){
                            continue;
                        }
                    }
                }
                countM[o1->start[j]-start]++;
            }
            for(j=0; j<o2->l; j++) {
                if(strand!=2){
                    if(ifp2->hdr->version & BM_STRAND){
                        if(strand != o2->strand[j]){
                            continue;
                        }
                    }
                }
                countM[o2->start[j]-start]++;
            }
            for(j=0; j<o1->l; j++) {
                if(countM[o1->start[j]-start] == 2){
                    for(k=lociK; k<o2->l; k++) {
                        if(o1->start[j] != o2->start[k]){
                            continue;
                        }
                        lociK = k;
                        printf("%s\t%ld", chrom, o1->start[j]);
                        if(ifp1->hdr->version & BM_CONTEXT)
                            printf("\t%s", context_str[o1->context[j]]);
                        if(ifp1->hdr->version & BM_STRAND)
                            printf("\t%s", strand_str[o1->strand[j]]);

                        printf("\t%f", o1->value[j]);
                        if(ifp1->hdr->version & BM_COVER)
                            printf("\t%"PRIu16"", o1->coverage[j]);
                        
                        printf("\t%f", o2->value[k]);
                        if(ifp2->hdr->version & BM_COVER)
                            printf("\t%"PRIu16"", o2->coverage[k]);

                        if(ifp1->hdr->version & BM_ID)
                            printf("\t%s", o1->entryid[j]);
                        printf("\n");
                    }
                }
            }

        }
        bmDestroyOverlappingIntervals(o1);
        bmDestroyOverlappingIntervals(o2);
    }

    //free(chrom);
    free(countM);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    return 1;
}

int bm_overlap_mul(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand){
    fprintf(stderr, "process region %s %d %d\n", chrom, start, end);
    int slen = 1, i =0, j = 0;
    int* countM = malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    
    //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
    if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
    int total = 0, loci = 0;
    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) goto error;
        //fprintf(stderr, "--- version --- %d\n", ifp1->hdr->version);
        if(o1->l) {
            for(j=0; j<o1->l; j++) {
                loci = o1->start[j]-start;
                if(strand!=2){
                    if(ifp1s[i]->hdr->version & BM_STRAND){
                        if(strand != o1->strand[j]){
                            continue;
                        }
                    }
                }
                countM[loci]++;
                if(countM[loci] == sizeifp){
                    total++;
                }
            }
        }
        bmDestroyOverlappingIntervals(o1);
    }
    if(total == 0){
        free(countM);
        return 0;
    }

    char **printmr = malloc(sizeof(char*)*(end-start+1));
    for(i=0;i<end-start+1;i++){
        printmr[i] = malloc(200+sizeifp*20);
    }
    char *tempchar = malloc(20);

    int printed = 0;
    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) goto error;
        if(o1->l){
            for(j=0; j<o1->l; j++) {
                loci = o1->start[j]-start;
                if(countM[loci] == sizeifp){
                    if(i == 0){
                        sprintf(printmr[loci],"%s\t%ld", chrom, o1->start[j]);
                        if(ifp1s[i]->hdr->version & BM_CONTEXT){
                            sprintf(tempchar,"\t%s", context_str[o1->context[j]]);
                            strcat(printmr[loci], tempchar);
                        }
                        if(ifp1s[i]->hdr->version & BM_STRAND){
                            sprintf(tempchar,"\t%s", strand_str[o1->strand[j]]);
                            strcat(printmr[loci], tempchar);
                        }

                    }

                    if(ifp1s[i]->hdr->version & BM_COVER){
                        sprintf(tempchar,"\t%d", (int)((double)o1->value[j]*o1->coverage[j] + 0.5));
                        strcat(printmr[loci], tempchar);
                        sprintf(tempchar,"\t%"PRIu16"", o1->coverage[j]);
                        strcat(printmr[loci], tempchar);
                    }else{
                        sprintf(tempchar,"\t%f", o1->value[j]);
                        strcat(printmr[loci], tempchar);
                    }
                    
                    if(i==sizeifp-1){
                        if(ifp1s[i]->hdr->version & BM_ID){
                            sprintf(tempchar,"\t%s", o1->entryid[j]);
                            strcat(printmr[loci], tempchar);
                        }
                    }
                }
            }

        }
        bmDestroyOverlappingIntervals(o1);
    }

    for(i=0;i<end-start+1;i++){
        if(countM[i] == sizeifp){
            printf("%s\n",printmr[i]);
        }
    }

    //free(chrom);
    for(i=0;i<end-start+1;i++){
        free(printmr[i]);
    }
    free(tempchar);
    free(printmr);
    free(countM);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    return 1;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE *Handle;
    Handle = fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		fprintf(stderr, "File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

void bmPrintHdr(binaMethFile_t *bm) {
    uint64_t i;
    int64_t i64;
    fprintf(stderr, "Ver code:    %"PRIu16"\n", bm->hdr->version);
    if(bm->hdr->version & BM_END) fprintf(stderr, "BM_END:    yes\n");
    else fprintf(stderr, "BM_END:    no\n");
    if(bm->hdr->version & BM_COVER) fprintf(stderr, "BM_COVER:    yes\n");
    else fprintf(stderr, "BM_COVER:    no\n");
    if(bm->hdr->version & BM_CONTEXT) fprintf(stderr, "BM_CONTEXT:    yes\n");
    else fprintf(stderr, "BM_CONTEXT:    no\n");
    if(bm->hdr->version & BM_STRAND) fprintf(stderr, "BM_STRAND:    yes\n");
    else fprintf(stderr, "BM_STRAND:    no\n");
    if(bm->hdr->version & BM_ID) fprintf(stderr, "BM_ID:    yes\n");
    else fprintf(stderr, "BM_ID:    no\n");
    fprintf(stderr, "Levels:     %"PRIu16"\n", bm->hdr->nLevels);
    //fprintf(stderr, "ctOffset:   0x%"PRIx64"\n", bm->hdr->ctOffset);
    //fprintf(stderr, "dataOffset: 0x%"PRIx64"\n", bm->hdr->dataOffset);
    fprintf(stderr, "indexOffset:        0x%"PRIx64"\n", bm->hdr->indexOffset);
    //fprintf(stderr, "sqlOffset:  0x%"PRIx64"\n", bm->hdr->sqlOffset);
    //fprintf(stderr, "summaryOffset:      0x%"PRIx64"\n", bm->hdr->summaryOffset);
    fprintf(stderr, "bufSize:    %"PRIu32"\n", bm->hdr->bufSize);
    fprintf(stderr, "extensionOffset:    0x%"PRIx64"\n", bm->hdr->extensionOffset);

    if(bm->hdr->nLevels) {
        fprintf(stderr, "	i	level	data	index\n");
    }
    for(i=0; i<bm->hdr->nLevels; i++) {
        fprintf(stderr, "\t%"PRIu64"\t%"PRIu32"\t%"PRIx64"\t%"PRIx64"\n", i, bm->hdr->zoomHdrs->level[i], bm->hdr->zoomHdrs->dataOffset[i], bm->hdr->zoomHdrs->indexOffset[i]);
    }

    fprintf(stderr, "nBasesCovered:      %"PRIu64"\n", bm->hdr->nBasesCovered);
    fprintf(stderr, "sumData:      %f\n", bm->hdr->sumData);
    fprintf(stderr, "minVal:     %f\n", bm->hdr->minVal);
    fprintf(stderr, "maxVal:     %f\n", bm->hdr->maxVal);
    //fprintf(stderr, "sumData:    %f\n", bm->hdr->sumData);
    //fprintf(stderr, "sumSquared: %f\n", bm->hdr->sumSquared);

    //Chromosome idx/name/length
    if(bm->cl) {
        fprintf(stderr, "Chromosome List\n");
        fprintf(stderr, "  idx\tChrom\tLength (bases)\n");
        for(i64=0; i64<bm->cl->nKeys; i64++) {
            fprintf(stderr, "  %"PRIu64"\t%s\t%"PRIu32"\n", i64, bm->cl->chrom[i64], bm->cl->len[i64]);
        }
    }
}

void bmPrintIndexNode(bmRTreeNode_t *node, int level) {
    uint16_t i;
    if(!node) return;
    for(i=0; i<node->nChildren; i++) {
        if(node->isLeaf) {
            printf("  %i\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t0x%"PRIx64"\t%"PRIu64"\n", level,\
                node->chrIdxStart[i], \
                node->baseStart[i], \
                node->chrIdxEnd[i], \
                node->baseEnd[i], \
                node->dataOffset[i], \
                node->x.size[i]);
        } else {
            printf("  %i\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t0x%"PRIx64"\tNA\n", level,\
                node->chrIdxStart[i], \
                node->baseStart[i], \
                node->chrIdxEnd[i], \
                node->baseEnd[i], \
                node->dataOffset[i]);
            bmPrintIndexNode(node->x.child[i], level+1);
        }
    }
}

void bmPrintIndexTree(binaMethFile_t *fp) {
    printf("\nIndex tree:\n");
    printf("nItems:\t%"PRIu64"\n", fp->idx->nItems);
    printf("chrIdxStart:\t%"PRIu32"\n", fp->idx->chrIdxStart);
    printf("baseStart:\t%"PRIu32"\n", fp->idx->baseStart);
    printf("chrIdxEnd:\t%"PRIu32"\n", fp->idx->chrIdxEnd);
    printf("baseEnd:\t%"PRIu32"\n", fp->idx->baseEnd);
    printf("idxSize:\t%"PRIu64"\n", fp->idx->idxSize);
    printf("  level\tchrIdxStart\tbaseStart\tchrIdxEnd\tbaseEnd\tchild\tsize\n");
    bmPrintIndexNode(fp->idx->root, 0);
}

char *fastStrcat(char *s, char *t)
{
	assert(s != NULL && t != NULL);

	while (*s != '\0')
		s++;

	//while ((*s++ = *t++) != '\0');
	while(1){
		if(*t == '\0'){ //  || *t == '\n'  || *t == '\r'
			*s++ = '\0';
			break;
		}
		*s++ = *t++;
	}

	return --s;
}

void onlyexecuteCMD(const char *cmd, const char *errorinfor)
{
    char ps[1024]={0};
    FILE *ptr;
    strcpy(ps, cmd);
    //fprintf(stderr, "[MM] %s\n", cmd);
    ptr=popen(ps, "w");

    if(ptr==NULL)
    {
        fprintf(stderr, "\n%s\n", errorinfor);
        exit(0);
    }
    pclose(ptr);
    ptr = NULL;
}

size_t get_executable_path( char* processdir,char* processname, size_t len)
{
	char* path_end;
	if(readlink("/proc/self/exe", processdir,len) <=0)
		return -1;
	path_end = strrchr(processdir, '/');
	if(path_end == NULL)
		return -1;
	++path_end;
	strcpy(processname, path_end);
	*path_end = '\0';
	return (size_t)(path_end - processdir);
}

unsigned long get_chr_len(binaMethFile_t *bm, char* chrom){
    int64_t i64;
    unsigned long chrlen = 0;
    //Chromosome idx/name/length
    if(bm->cl) {
        for(i64=0; i64<bm->cl->nKeys; i64++) {
            //fprintf(stderr, "  %"PRIu64"\t%s\t%"PRIu32"\n", i64, bm->cl->chrom[i64], bm->cl->len[i64]);
            if(strcmp(bm->cl->chrom[i64], chrom) == 0) {
                chrlen = bm->cl->len[i64];
                fprintf(stderr, "Chromosome List %s %ld\n", chrom, chrlen);
            }
        }
    }
    return chrlen;
}
