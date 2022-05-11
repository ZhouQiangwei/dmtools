/*
* This tools is used for make bigmeth file and process bigmeth file.
* main function contains:
* bam2bm calculate methleval in loci and region with bam file, output bm file.
* view   view bigmeth as txt format, 20211101
* mr2mbw  convert methratio file from batmeth2 to bigmeth, 20211103
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
* ./exampleWrite mr2mbw -g ~/practice/Genome/hg38/hg38.chr.fa.len -E -C -m test.f -S --Cx -o test.bm -r chr1:0-100,chr1:16766-16890
* gcc test/exampleWrite.c -o exampleWrite -I. -L. -lBigWig -Wl,-rpath /public/home/qwzhou/software_devp/batmeth2-bwa/src/bmtools/
*/
#include "bigWig.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

FILE* File_Open(const char* File_Name,const char* Mode);
char *strand_str[] = {"+", "-", "."};
char *context_str[] = {"C", "CG", "CHG", "CHH"};
int main_view_all(bigWigFile_t *fp, FILE* outfileF, char *outformat, bigWigFile_t *ofp, char* filterchrom);
int main_view(bigWigFile_t *ifp, char *region, FILE* outfileF, char *outformat, bigWigFile_t *ofp);
int main_view_bedfile(char *inbmF, char *bedfile, int type, FILE* outfileF, char *outformat, bigWigFile_t *ofp);
int calchromstats(char *inbmfile, char *method, int chromstep, int stepoverlap, uint8_t strand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calregionstats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calregionstats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calbodystats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int calbodystats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh);
int bw_overlap_region(char *inbmF1, char *inbmF2, char *region, uint8_t pstrand);
int bw_overlap_region_mul(char *inbmFs, char *region, uint8_t pstrand);
int bw_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand);
int bw_overlap_file_mul(char *inbmF1, char *bedfile);
int bw_overlap_file(char *inbmF1, char *inbmF2, char *bedfile);
int bw_overlap_all_mul(char *inbmFs, uint8_t pstrand);
int bw_overlap_mul(bigWigFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand);
int bw_overlap(bigWigFile_t *ifp1, bigWigFile_t *ifp2, char *chrom, int start, int end, uint8_t strand);
int calprofile(char *inbmfile, int upstream, int downstream, double profilestep, double bodyprofilestep, double bodyprofilemovestep, double profilemovestep, char *bedfile, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX);
int calprofile_gtf(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *gtffile, int format, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX);
void getcontext(uint8_t context, char* str_context);
void delete_char2(char* str,char target,char target2);
void calbody_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF);
void calregion_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* strand, FILE* outfileF);
void calregion_weighted_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, uint16_t *countC, uint16_t *countCT);
int main_view_file(bigWigFile_t *ifp, char *bedfile, FILE* outfileF, char *outformat, bigWigFile_t *ofp);
void mbwfileinit(bigWigFile_t *ofp, bigWigFile_t *ifp, char* outfile, int zoomlevel);
void bwPrintHdr(bigWigFile_t *bw);
void bwPrintIndexNode(bwRTreeNode_t *node, int level);

#define MAX_LINE_PRINT 1000000
#define MAX_BUFF_PRINT 20000000
const char* Help_String_main="Command Format :  bmtools <mode> [opnions]\n"
		"\nUsage:\n"
        "\t  [mode]         mr2mbw view viewheader overlap regionstats bodystats profile chromstats\n\n"
        "\t  mr2mbw         convert txt meth file to mbw format\n"
        "\t  view           mbw format to txt meth\n"
        "\t  viewheader     view header of mbw file\n"
        "\t  overlap        overlap cytosine site with more than two mbw files\n"
        "\t  regionstats    calculate DNA methylation level of per region\n"
        "\t  bodystats      calculate DNA methylation level of body, upstream and downstream.\n"
        "\t  profile        calculate DNA methylation profile\n"
        "\t  chromstats     calculate DNA methylation level across chromosome\n";

const char* Help_String_mr2mbw="Command Format :  bmtools mr2mbw [opnions] -g genome.fa.fai -m methratio.txt -o outmeth.mbw\n"
		"\nUsage: bmtools mr2mbw -C -S --Cx -E -g genome.fa.fai -m meth.txt -o meth.mbw\n"
        "\t [mr2mbw] mode paramaters, required\n"
		"\t-g                    chromosome size file.\n"
        "\t-m                    methratio file\n"
        "\t-o|--outmbw           output mbigwig file\n"
        "\t [mr2mbw] mode paramaters, options\n"
        "\t-C                    print coverage\n"
        "\t-S                    print strand\n"
        "\t--Cx                  print context\n"
        "\t-E                    print end\n"
        "\t--Id                  print ID\n"
        "\t--CF                  coverage filter, >=[int], default 4.\n"
        "\t--sort Y/N            make chromsize file and meth file in same coordinate, default Y\n"
        "\t--zl                  The maximum number of zoom levels. [1-10]\n"
        "\t-f                    file format. methratio, bedmethyl or bedsimple [default methratio]\n"
        "\t  methratio           chrom start strand context meth_reads cover\n"
        "\t  bedmethyl           chrom start end name * strand * * * coverage meth_reads\n"
        "\t  bedsimple           chrom start end id strand context meth_reads coverage\n"
        "\t--pcontext            CG/CHG/CHH/C, needed when bedmethyl format, default C\n"
        "\t--context             CG/CHG/CHH/ALL, only convert provide context in methratio file or bedsimple, default ALL\n"
        "\tNote. meth ratio file must be sorted by chrom and coordinate. ex. sort -k1,1 -k2,2n\n"
		"\t-h|--help";

const char* Help_String_view="Command Format :  bmtools view [opnions] -i meth.mbw\n"
		"\nUsage:\n"
        "\t [view] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t [view] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by semicolon. chr1:1-2900;chr2:1-200;\n"
        "\t--bed                 bed file for view, format: chrom start end (strand).\n"
        "\t--strand              [0/1/2] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information\n"
        "\t--context             [0/1/2/3] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context.\n"
        "\t--mincover            >= minumum coverage show, default: 0\n"
        "\t--maxcover            <= maximum coverage show, default: 10000\n"
		"\t-h|--help";

const char* Help_String_viewheader="Command Format :  bmtools viewheader -i meth.mbw\n"
		"\nUsage:\n"
        "\t [view] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
		"\t-h|--help";

const char* Help_String_overlap="Command Format :  bmtools overlap [opnions] -i meth1.mbw -i2 meth2.mbw\n"
		"\nUsage:\n"
        "\t [overlap] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t-i2                   input mbigwig file2\n"
        "\t [overlap] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by semicolon. chr1:1-2900;chr2:1-200;\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t [overlap] mode paramaters\n"
        "\t--bmfiles             input mbigwig files, seperated by comma. This is no need if you provide -i and -i2.\n"
		"\t-h|--help";

const char* Help_String_regionstats="Command Format :  bmtools regionstats [opnions] -i mbw --bed bedfile\n"
		"\nUsage:\n"
        "\t [regionstats] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
        "\t [regionstats] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by semicolon. chr1:1-2900;chr2:1-200,+;\n"
        "\t--method              weighted/ mean\n"
        "\t--strand              [0/1/2/3] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information, 3 calculate and print strand meth level seperately\n"
        "\t--context             [0/1/2/3/4] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context, 4 calculate and print strand meth level seperately\n"
        "\t--printcoverage       [0/1] print countC and coverage instead of methratio. [0]\n"
        "\t--print2one           [int] print all the countC and coverage results of C/CG/CHG/CHH context methylation to same file, only valid when --printcoverage 1. 0 for no, 1 for yes. [0]\n"
//        "\t--mincover            >= minumum coverage show, default: 0\n"
//        "\t--maxcover            <= maximum coverage show, default: 10000\n"
		"\t-h|--help";

const char* Help_String_bodystats="Command Format :  bmtools bodystats [opnions] -i mbw --bed bedfile\n"
		"\nUsage:\n"
        "\t [bodystats] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
        "\t [bodystats] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t-r                    region for view, can be seperated by semicolon. chr1:1-2900;chr2:1-200,+;\n"
        "\t--method              weighted/ mean\n"
        "\t--regionextend        also calculate DNA methylation level of upstream and downstream N-bp window. default 2000.\n"
        "\t--strand              [0/1/2/3] strand for show, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' all information, 3 calculate and print strand meth level seperately\n"
        "\t--context             [0/1/2/3/4] context for show, 0 represent 'C/ALL' context, 1 'CG' context, 2 'CHG' context, 3 'CHH' context, 4 calculate and print strand meth level seperately\n"
		"\t--printcoverage       [0/1] also print countC and coverage of body instead of methratio. [0]\n"
        "\t--print2one           [int] print all the countC and coverage results of C/CG/CHG/CHH context methylation to same file, only valid when --printcoverage 1. 0 for no, 1 for yes. [0]\n"
        "\t-h|--help";

const char* Help_String_profile="Command Format :  bmtools profile [opnions] -i meth.mbw --bed bedfile\n"
		"\nUsage:\n"
        "\t [profile] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
		"\t [profile] mode paramaters, options\n"
        "\t-o                    output file [stdout]\n"
        "\t--regionextend        region extend for upstream and downstram, [2000]\n"
        "\t--profilestep         [double] step mean bin size for chromosome region, default: 0.02 (2%)\n"
        "\t--profilemovestep     [double] step move, default: 0.01, if no overlap, please define same as --profilestep\n"
        "\t--profilemode         calculate the methylation matrix mode of every region or gene. 0 for gene and flanks mode, 1 for tss and flanks, 2 for tts and flanks, 3 for region center and flanks. [0]\n"
        "\t--bodyX               [double] the gene body bin size is N multiple of the bin size of the upstream and downstream extension region. [1]\n"
        "\t--matrixX             [int] the bin size is N times of the bin size of profile, so the bin size should be N*profilestep [5], please note N*profilestep must < 1 and N must >= 1\n"
        "\t--print2one           [int] print all the matrix results of C/CG/CHG/CHH context methylation to same file. 0 for no, 1 for yes. [0]\n"
        "\t-h|--help";

const char* Help_String_chromstats="Command Format :  bmtools chromstats [opnions] -i meth.mbw\n"
		"\nUsage:\n"
        "\t [chromstats] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
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

int main(int argc, char *argv[]) {
    bigWigFile_t *fp = NULL;
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

    char *mrformat = malloc(100);
    strcpy(mrformat, "methratio");
    char *pcontext = malloc(10);
    strcpy(pcontext, "C");
    char *filtercontext = malloc(10);
    strcpy(filtercontext, "ALL");
    // for gene file
    int upstream = 2000, downstream = 2000;
    double profilestep = 0.02, profilemovestep = 0.01;
    double bodyprofilestep = 0.02, bodyprofilemovestep = 0.01;
    unsigned int mcover_cutoff = 4; unsigned long TotalC = 0;
    int zoomlevel = 5;
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
        //mr2mbw/ view/ overlap/ regionstats/ profile/ chromstats
        if(strcmp(mode, "mr2mbw") == 0){
           fprintf(stderr, "%s\n", Help_String_mr2mbw); 
        }else if(strcmp(mode, "view") == 0){
           fprintf(stderr, "%s\n", Help_String_view); 
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
        }else{
            fprintf(stderr, "%s\n", Help_String_main);
        }
        exit(0);
    }else{
        fprintf(stderr, "Please define mode!!!\n");
        fprintf(stderr, "%s\n", Help_String_main);
        exit(0);
    }
    for(i=0; i< argc; i++){
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
        }else if(strcmp(argv[i], "--outmbw") == 0){
            outbmfile = malloc(100);
            strcpy(outbmfile, argv[++i]);
        }else if(strcmp(argv[i], "--outformat") == 0){
            strcpy(outformat, argv[++i]);
        }else if(strcmp(argv[i], "-i") == 0){
            strcpy(inbmfile, argv[++i]);
        }else if(strcmp(argv[i], "--CF") == 0){
            mcover_cutoff = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--bmfiles") == 0){
            strcpy(inbmfiles, argv[++i]);
            inbm_mul = 1;
        }else if(strcmp(argv[i], "-i2") == 0){
            strcpy(bmfile2, argv[++i]);
        }else if(strcmp(argv[i], "-r") == 0){
            region = malloc(1000);
            strcpy(region, argv[++i]);
        }else if(strcmp(argv[i], "--method") == 0){
            strcpy(method, argv[++i]);
        }else if(strcmp(argv[i], "--chromstep") == 0){
            chromstep = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--chrom") == 0){
            filterchrom = malloc(200);
            strcpy(filterchrom, argv[++i]);
        }else if(strcmp(argv[i], "--stepmove") == 0){
            stepoverlap = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--profilestep") == 0){
            profilestep = atof(argv[++i]);
            assert(profilestep>0 && profilestep<1);
        }else if(strcmp(argv[i], "--regionextend") == 0){
            regionextend = atoi(argv[++i]);
            assert(regionextend>0);
        }else if(strcmp(argv[i], "--profilemovestep") == 0){
            profilemovestep = atof(argv[++i]);
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
            gfffile = malloc(100);
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
        }else if(strcmp(argv[i], "--zl") == 0){
            zoomlevel = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--fcontext") == 0){
            strcpy(filtercontext, argv[++i]);
        }else if(strcmp(argv[i], "--regionextend") == 0){
            regionextend = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--mincover") == 0){
            mincover = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--maxcover") == 0){
            maxcover = atoi(argv[++i]);
        }else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--out") == 0){
            if(strcmp(mode, "mr2mbw") == 0){
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

    if(strcmp(mode, "mr2mbw") == 0){
        fprintf(stderr, "mr file format %s\n", mrformat);
        if(chromlenf_yes==0){
            fprintf(stderr, "please provide chrome size file with -g paramater\n");
            exit(0);
        }
        FILE *methF = File_Open(methfile, "r"); 
        char *chrom = malloc(50); char *old_chrom = malloc(50);
        //int MAX_CHROM = 10000;
        char *PerLine = malloc(200); 
        //char **chromsArray = malloc(sizeof(char*)*MAX_CHROM);
        unsigned long chrprintL = 0;
        if(strcmp(sortY, "Y") == 0){
            fprintf(stderr, "obtained chromosome order in meth ratio file ... \n");
            while(fgets(PerLine,200,methF)!=0){
                if(PerLine[0] == '#') continue; // remove header #
                //fprintf(stderr, "%s\n", PerLine);
                if(strcmp(mrformat, "methratio") == 0 || strcmp(mrformat, "bedmethyl") == 0 || strcmp(mrformat, "bedsimple") == 0){
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
        while(fgets(PerLine,200,ChromF)!=NULL){
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
        if(bwInit(1<<17) != 0) {
            fprintf(stderr, "Received an error in bwInit\n");
            return 1;
        }

        fp = bwOpen(outbmfile, NULL, "w");
        fp->type = write_type;
        if(!fp) {
            fprintf(stderr, "An error occurred while opening example_output.bw for writingn\n");
            return 1;
        }
        
        //Allow up to 10 zoom levels, though fewer will be used in practice
        if(bwCreateHdr(fp, zoomlevel)) {
            fprintf(stderr, "== bwCreateHdr ==\n");
            goto error;
        }

        //Create the chromosome lists
        fp->cl = bwCreateChromList(chroms, chrLens, printL); //2
        if(!fp->cl) {
            fprintf(stderr, "== bwCreateChromList ==\n");
            goto error;
        }

        //Write the header
        if(bwWriteHdr(fp)) {
            fprintf(stderr, "== bwWriteHdr ==\n");
            goto error;
        }
        
        //Some example methlevel
        if(DEBUG>1) fprintf(stderr, "====HHH type %d\n", fp->type);
        methF = File_Open(methfile, "r");
        strcpy(old_chrom, "NN"); 
        char *strand = malloc(2), *context = malloc(10), *nameid = malloc(100);
        unsigned start=0, end = 0; unsigned int coverC,coverage=0; float value=0;
        printL = 0;
        strcpy(context, pcontext);
        char *decide = malloc(10);
        while(fgets(PerLine,200,methF)!=0){
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
                //bedmethyl2mbw
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
            }else{
                fprintf(stderr, "Unexpected mr file format!!!\n");
                exit(0);
            }
            if(coverage<mcover_cutoff) continue;
            if(strcmp(filtercontext, "ALL") != 0 && strcmp(filtercontext, context) != 0){
                continue;
            }

            if(strcmp(old_chrom, chrom)!=0){
                if(printL > 1){
                    if(bwAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1))  {
                        fprintf(stderr, "bwAppendIntervals -1\n");
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
                if(DEBUG>-1) fprintf(stderr,"%d XXX %s %d %d %f %d %d %d\n", printL, chromsUse[printL], starts[printL], ends[printL], values[printL], coverages[printL], strands[printL], context[printL]);
                int response = bwAddIntervals(fp, chromsUse, starts, ends, values, coverages, strands, contexts, 
                entryid, 1);
                fprintf(stderr, "Processing %s chromosome.\n", chrom);
                if(response) {
                    fprintf(stderr, "bwAddIntervals 0\n");
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
                if(bwAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) {
                    fprintf(stderr, "bwAppendIntervals 2\n");
                    goto error;
                }
                TotalC+=printL;
                printL = 0;
            }
        } // end read me file
        if(printL > 1) {
            if(DEBUG>1) fprintf(stderr, "--last print %d %d %d\n", starts[printL-1], ends[printL-1], printL-1);
            if(bwAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) {
                fprintf(stderr, "bwAppendIntervals 3\n");
                goto error;
            }
            TotalC+=printL;
            printL = 0; 
        }
        fprintf(stderr, "\nprocess %ld cytosine site\n", TotalC);
        //Add a new block of entries with a span. Since bwAdd/AppendIntervals was just used we MUST create a new block
    //    if(bwAddIntervalSpans(fp, "1", starts+6, 20, values+6, 3)) goto error;
        //We can continue appending similarly formatted entries
    //    if(bwAppendIntervalSpans(fp, starts+9, values+9, 3)) goto error;

        //Add a new block of fixed-step entries
    //    if(bwAddIntervalSpanSteps(fp, "1", 900, 20, 30, values+12, 3)) goto error;
        //The start is then 760, since that's where the previous step ended
    //    if(bwAppendIntervalSpanSteps(fp, values+15, 3)) goto error;

        //Add a new chromosome
    //    chromsUse[0] = "2";
    //    chromsUse[1] = "2";
    //    chromsUse[2] = "2";
    //    if(bwAddIntervals(fp, chromsUse, starts, ends, values, 3)) goto error;

        //Closing the file causes the zoom levels to be created
        //if(DEBUG>0) 
        if(DEBUG>0) fprintf(stderr, "bm close1111 ----- \n");
        bwClose(fp);
        if(DEBUG>0) fprintf(stderr, "bm close22222 ---===--- \n");
        bwCleanup();
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
    if(strcmp(mode, "view")==0){
        if(DEBUG>0) fprintf(stderr, "bm view\n");
        //char region[] = "chr1:0-100,chr1:16766-16830";
        uint32_t type = BMtype(inbmfile, NULL);
        bigWigFile_t *ifp = NULL;
        ifp = bwOpen(inbmfile, NULL, "r");
        ifp->type = ifp->hdr->version;
        FILE *outfp_mbw = NULL;
        bigWigFile_t *ofp = NULL;
        if(strcmp(outformat, "txt") == 0){
            if(outfile) {
                outfp_mbw = File_Open(outfile,"w");
            }else {
                outfp_mbw = stdout;
            }
        }else if(strcmp(outformat, "mbw") == 0){
            fprintf(stderr, "output mbw file\n");
            //mbwfileinit(ofp, ifp, outfile, zoomlevel); // why not valid????
            if(bwInit(1<<17) != 0) {
                fprintf(stderr, "Received an error in bwInit\n");
                return;
            }
            ofp = bwOpen(outfile, NULL, "w");
            ofp->type = ifp->type;
            if(!ofp) {
                fprintf(stderr, "An error occurred while opening mbw for writingn\n");
                return;
            }
            //Allow up to 10 zoom levels, though fewer will be used in practice
            if(bwCreateHdr(ofp, zoomlevel)) {
                fprintf(stderr, "== bwCreateHdr ==\n");
                return;
            }
            //Create the chromosome lists
            ofp->cl = bwCreateChromList_ifp(ifp); //2
            if(!ofp->cl) {
                fprintf(stderr, "== bwCreateChromList ==\n");
                return;
            }
            //Write the header
            if(bwWriteHdr(ofp)) {
                fprintf(stderr, "== bwWriteHdr ==\n");
                return;
            }
        }
        
        if(!region && !bedfile){
            fprintf(stderr, "[Mode] ------------- View all meth\n");
            main_view_all(ifp, outfp_mbw, outformat, ofp, filterchrom);
        }else if(region){
            fprintf(stderr, "[Mode] ------------- View region %s\n", region);
            main_view(ifp, region, outfp_mbw, outformat, ofp);
            free(region);
        }else if(bedfile){
            fprintf(stderr, "[Mode] ------------- View bedfile %s\n", bedfile);
            main_view_file(ifp, bedfile, outfp_mbw, outformat, ofp);
            free(bedfile);
        }else{
            fprintf(stderr, "\nplease provide -r or --bed!!!\n");
        }
        fprintf(stderr, "Done and free mem\n");
        free(inbmfile);
        bwClose(ifp);
        if(outfile) free(outfile);
        if(filterchrom) free(filterchrom);
        if(strcmp(outformat, "txt") == 0){
            fclose(outfp_mbw);
        }else if(strcmp(outformat, "mbw") == 0){
            bwClose(ofp);
            bwCleanup();
        }
        return 0;
    }

    if(strcmp(mode, "viewheader")==0){
        if(DEBUG>0) fprintf(stderr, "bm viewheader\n");
        uint32_t type = BMtype(inbmfile, NULL);
        bigWigFile_t *ifp = NULL;
        ifp = bwOpen(inbmfile, NULL, "r");
        ifp->type = ifp->hdr->version;
        bwPrintHdr(ifp);
        //bwPrintIndexTree(ifp);
        free(inbmfile);
        bwClose(ifp);
        return 0;
    }

    //overlap
    if(strcmp(mode, "overlap")==0){
        fprintf(stderr, "[Mode] ------------- overlap %s %s\n", inbmfile, bmfile2);
        if(inbm_mul == 1){
            if(!region && !bedfile){
                bw_overlap_all_mul(inbmfiles, filter_strand);
            }else if(region){
                bw_overlap_region_mul(inbmfiles, region, filter_strand);
                free(region);
            }else if(bedfile){
                bw_overlap_file_mul(inbmfiles, bedfile);
                free(bedfile);
            }else{
                fprintf(stderr, "\nplease provide -r or --bed!!!\n");
            }
            free(inbmfiles); 
        }else{
            if(!region && !bedfile){
                bw_overlap_all(inbmfile, bmfile2, 1, 1, filter_strand);
            }else if(region){
                bw_overlap_region(inbmfile, bmfile2, region, filter_strand);
                free(region);
            }else if(bedfile){
                bw_overlap_file(inbmfile, bmfile2, bedfile);
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
            fclose(outfp);
            if(outfile && print2one == 0) { fclose(outfp_cg); fclose(outfp_chg); fclose(outfp_chh); }
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
        if(outfile) {
            char* outfile_aver_temp = malloc(sizeof(char)*200);
            strcpy(outfile_aver_temp, outfile); strcat(outfile_aver_temp, ".chrom");
            outfp_mean = File_Open(outfile_aver_temp,"w");
            free(outfile_aver_temp);
        }else {
            outfp_mean = stdout;
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
        fclose(outfp_mean); 
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

        //bodyX
        bodyprofilestep = profilestep * bodyX;
        bodyprofilemovestep = profilemovestep * bodyX;

        if(bedfile){
            calprofile(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, bedfile, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX);
            free(bedfile);
        }
        else if(gtffile){
            calprofile_gtf(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, gtffile, 1, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX);
            free(gtffile);
        }else if(gfffile){
            calprofile_gtf(inbmfile, upstream, downstream, profilestep, profilemovestep, bodyprofilestep, bodyprofilemovestep, gfffile, 2, filter_context, outfp, outfp_cg, outfp_chg, outfp_chh, outfp_aver, profilemode, matrixX);
            free(gfffile);
        }else{
            fprintf(stderr, "\nplease provide -r, --bed, --gtf or --gff!!!\n");
        }
        
        fclose(outfp);
        if(outfile && print2one == 0) { fclose(outfp_cg); fclose(outfp_chg); fclose(outfp_chh); }
        fclose(outfp_aver); 
        if(outfile) free(outfile);
        return 0;
    }

    if(outfile) free(outfile);
    return 1;
error:
    fprintf(stderr, "Received an error in process!\n");
    bwClose(fp);
    bwCleanup();
    return -1;
}

void mbwfileinit(bigWigFile_t *ofp, bigWigFile_t *ifp, char* outfile, int zoomlevel){
    if(bwInit(1<<17) != 0) {
        fprintf(stderr, "Received an error in bwInit\n");
        return;
    }
    ofp = bwOpen(outfile, NULL, "w");
    ofp->type = ifp->type;
    if(!ofp) {
        fprintf(stderr, "An error occurred while opening mbw for writingn\n");
        return;
    }
    //Allow up to 10 zoom levels, though fewer will be used in practice
    if(bwCreateHdr(ofp, zoomlevel)) {
        fprintf(stderr, "== bwCreateHdr ==\n");
        return;
    }

    //Create the chromosome lists
    ofp->cl = ifp->cl; //2
    if(!ofp->cl) {
        fprintf(stderr, "== bwCreateChromList ==\n");
        return;
    }

    //Write the header
    if(bwWriteHdr(ofp)) {
        fprintf(stderr, "== bwWriteHdr ==\n");
        return;
    }
}

double *Sregionstats(bigWigFile_t *fp, char *chrom, int start, int end, int splitN, uint32_t movestep, char *method, uint8_t strand, uint8_t context){
    double *stats = NULL;
    assert(splitN>0);
    //int i=0;
    if(strcmp(method, "mean")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, mean, strand, context);
    }else if(strcmp(method, "weighted")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, weighted, strand, context);
    }else if(strcmp(method, "dev")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, dev, strand, context);
    }else if(strcmp(method, "min")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, min, strand, context);
    }else if(strcmp(method, "max")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, max, strand, context);
    }else if(strcmp(method, "cover")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, cov, strand, context);
    }
    return stats;
}

//c cg chg chh
//idx0 c cg chg chh; idx1 c cg chg chh; ... etc
double *Sregionstats_array(bigWigFile_t *fp, char *chrom, int start, int end, int splitN, uint32_t movestep, char *method, uint8_t strand){
    double *stats = NULL;
    assert(splitN>0);
    //int i=0;
    if(strcmp(method, "mean")==0){
        stats = bwStats_array(fp, chrom, start, end, splitN, movestep, mean, strand);
    }else if(strcmp(method, "weighted")==0){
        stats = bwStats_array(fp, chrom, start, end, splitN, movestep, weighted, strand);
    }
    return stats;
}

//double *output = malloc(sizeof(double)*nBins*Tsize);
void *Sregionstats_array_count(bigWigFile_t *fp, char *chrom, int start, int end, int splitN, uint32_t movestep, char *method, uint8_t strand, uint16_t *countC, uint16_t *countCT){
    assert(splitN>0);
    int i=0, Tsize = 4;
    for(i=0;i<splitN*Tsize;i++){
        countC[i]=0;
        countCT[i]=0;
    }
    if(strcmp(method, "weighted")==0){
        bwStats_array_count(fp, chrom, start, end, splitN, movestep, weighted, strand, countC, countCT);
    }else {
        fprintf(stderr, "Unexpected method for weighted methylation count!");
        exit(0);
    }
    return;
}

int calchromstats(char *inbmfile, char *method, int chromstep, int stepoverlap, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    int i = 0, start = 0, end = chromstep; //, j = 0
    char* region = malloc(sizeof(char)*1000);
    int splitN = 1, Tsize = 4;
    uint16_t *countC = malloc(sizeof(uint16_t)*splitN*Tsize);
    uint16_t *countCT = malloc(sizeof(uint16_t)*splitN*Tsize);
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

    bwClose(fp);
    bwCleanup();
    free(countC); free(countCT);
    return 0;
}

void profile_print(bigWigFile_t *fp, char *chrom, int start, int end, int splitN, double profilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, FILE* outfileF){
    double *stats = Sregionstats(fp, chrom, start, end, splitN, (int)((end-start)*profilemovestep), "weighted", strand, context);
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

void profile_print_array(double *finalstats, int *finalcounts, bigWigFile_t *fp, char *chrom, int start, int end, int splitN, double profilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, FILE* outfileF){
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, (int)((end-start)*profilemovestep), "weighted", strand);
    int i,j,k; int Tsize = 4; int total_splitN = splitN*Tsize;
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
    char* print_context = malloc(sizeof(char)*5);
    if(stats) {
        if(strand == 1) {// -
            for(j=total_splitN-Tsize, k=0; j<total_splitN; j++,k++){
                if(!(context>=4 || k==context)){
                    continue;
                }
                getcontext(k, print_context);
                if(geneid[0]) fprintf(outfileF, "%s:%d-%d-%s\t%s\t%f", chrom, start, end, geneid, print_context, stats[j]);
                else fprintf(outfileF, "%s:%d-%d\t%s\t%f", chrom, start, end, print_context, stats[j]);
                //store for average mr
                if(!isnan(stats[j])) {
                    finalstats[k] += stats[j];
                    finalcounts[k]++;
                }
                if(splitN>1){
                    for(i=j-Tsize;i>=0;){
                        fprintf(outfileF, "\t%f", stats[i]);
                        if(!isnan(stats[i])) {
                            finalstats[total_splitN-i] += stats[i];
                            finalcounts[total_splitN-i]++;
                        }
                        i=i-Tsize;
                    }
                }
                fprintf(outfileF, "\n");
            }
        }else{
            //0 C
            for(j=0; j<Tsize; j++){
                if(!(context>=4 || j==context)){
                    continue;
                }
                getcontext(j, print_context);
                if(geneid[0]) fprintf(outfileF, "%s:%d-%d-%s\t%s\t%f", chrom, start, end, geneid, print_context, stats[j]);
                else fprintf(outfileF, "%s:%d-%d\t%s\t%f", chrom, start, end, print_context, stats[j]);
                if(!isnan(stats[j])) {
                    finalstats[j] += stats[j];
                    finalcounts[j]++;
                }
                if(splitN>1){
                    for(i=j+Tsize;i<total_splitN;){
                        fprintf(outfileF, "\t%f", stats[i]);
                        if(!isnan(stats[i])) {
                            finalstats[i] += stats[i];
                            finalcounts[i]++;
                        }
                        i+=Tsize;
                    }
                }
                fprintf(outfileF, "\n");
            }
        }
    }
    free(print_context);
}
uint32_t stored_buffer = 0;
void profile_print_array_buffer(double *finalstats, int *finalcounts, bigWigFile_t *fp, char *chrom, int processtart, int start, int end, int processend, int splitN, int bodysplitN, double profilemovestep, double bodyprofilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, char* printbuffer_c, char* printbuffer_cg, char* printbuffer_chg, char* printbuffer_chh, int matrixX){
    double *stats_up = Sregionstats_array(fp, chrom, processtart, start, splitN, (int)((start-processtart)*profilemovestep), "weighted", strand);
    double *stats_body = Sregionstats_array(fp, chrom, start, end, bodysplitN, (int)((end-start)*bodyprofilemovestep), "weighted", strand);
    double *stats_down = Sregionstats_array(fp, chrom, end, processend, splitN, (int)((processend-end)*profilemovestep), "weighted", strand);
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
                if(k==0) strcat(printbuffer_c, storetemp);
                else if(k==1) strcat(printbuffer_cg, storetemp);
                else if(k==2) strcat(printbuffer_chg, storetemp);
                else if(k==3) strcat(printbuffer_chh, storetemp);
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
                        if(k==0) strcat(printbuffer_c, storetemp);
                        else if(k==1) strcat(printbuffer_cg, storetemp);
                        else if(k==2) strcat(printbuffer_chg, storetemp);
                        else if(k==3) strcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
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
                    if(k==0) strcat(printbuffer_c, storetemp);
                    else if(k==1) strcat(printbuffer_cg, storetemp);
                    else if(k==2) strcat(printbuffer_chg, storetemp);
                    else if(k==3) strcat(printbuffer_chh, storetemp);

                }
                sprintf(storetemp, "\n");
                if(k==0) strcat(printbuffer_c, storetemp);
                else if(k==1) strcat(printbuffer_cg, storetemp);
                else if(k==2) strcat(printbuffer_chg, storetemp);
                else if(k==3) strcat(printbuffer_chh, storetemp);
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
                if(j==0) strcat(printbuffer_c, storetemp);
                else if(j==1) strcat(printbuffer_cg, storetemp);
                else if(j==2) strcat(printbuffer_chg, storetemp);
                else if(j==3) strcat(printbuffer_chh, storetemp);
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
                        if(j==0) strcat(printbuffer_c, storetemp);
                        else if(j==1) strcat(printbuffer_cg, storetemp);
                        else if(j==2) strcat(printbuffer_chg, storetemp);
                        else if(j==3) strcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
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
                    if(j==0) strcat(printbuffer_c, storetemp);
                    else if(j==1) strcat(printbuffer_cg, storetemp);
                    else if(j==2) strcat(printbuffer_chg, storetemp);
                    else if(j==3) strcat(printbuffer_chh, storetemp);

                }
                sprintf(storetemp, "\n");
                if(j==0) strcat(printbuffer_c, storetemp);
                else if(j==1) strcat(printbuffer_cg, storetemp);
                else if(j==2) strcat(printbuffer_chg, storetemp);
                else if(j==3) strcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
}

void profile_print_array_buffer1(double *finalstats, int *finalcounts, bigWigFile_t *fp, char *chrom, int start, int end, int splitN, double profilemovestep, char *method, uint8_t strand, char* geneid, uint8_t context, char* printbuffer_c, char* printbuffer_cg, char* printbuffer_chg, char* printbuffer_chh, int matrixX){
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, (int)((end-start)*profilemovestep), "weighted", strand);
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
                if(k==0) strcat(printbuffer_c, storetemp);
                else if(k==1) strcat(printbuffer_cg, storetemp);
                else if(k==2) strcat(printbuffer_chg, storetemp);
                else if(k==3) strcat(printbuffer_chh, storetemp);
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
                        if(k==0) strcat(printbuffer_c, storetemp);
                        else if(k==1) strcat(printbuffer_cg, storetemp);
                        else if(k==2) strcat(printbuffer_chg, storetemp);
                        else if(k==3) strcat(printbuffer_chh, storetemp);

                        i=i-Tsize*matrixX;
                    }
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
                    if(k==0) strcat(printbuffer_c, storetemp);
                    else if(k==1) strcat(printbuffer_cg, storetemp);
                    else if(k==2) strcat(printbuffer_chg, storetemp);
                    else if(k==3) strcat(printbuffer_chh, storetemp);

                }
                sprintf(storetemp, "\n");
                if(k==0) strcat(printbuffer_c, storetemp);
                else if(k==1) strcat(printbuffer_cg, storetemp);
                else if(k==2) strcat(printbuffer_chg, storetemp);
                else if(k==3) strcat(printbuffer_chh, storetemp);
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
                if(j==0) strcat(printbuffer_c, storetemp);
                else if(j==1) strcat(printbuffer_cg, storetemp);
                else if(j==2) strcat(printbuffer_chg, storetemp);
                else if(j==3) strcat(printbuffer_chh, storetemp);
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
                        if(j==0) strcat(printbuffer_c, storetemp);
                        else if(j==1) strcat(printbuffer_cg, storetemp);
                        else if(j==2) strcat(printbuffer_chg, storetemp);
                        else if(j==3) strcat(printbuffer_chh, storetemp);

                        i+=Tsize*matrixX;
                    }
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
                    if(j==0) strcat(printbuffer_c, storetemp);
                    else if(j==1) strcat(printbuffer_cg, storetemp);
                    else if(j==2) strcat(printbuffer_chg, storetemp);
                    else if(j==3) strcat(printbuffer_chh, storetemp);

                }
                sprintf(storetemp, "\n");
                if(j==0) strcat(printbuffer_c, storetemp);
                else if(j==1) strcat(printbuffer_cg, storetemp);
                else if(j==2) strcat(printbuffer_chg, storetemp);
                else if(j==3) strcat(printbuffer_chh, storetemp);
            }
        }
    }
    //fprintf(stderr, "%s %d\n", printbuffer, strlen(printbuffer));
    free(print_context);
    free(storetemp);
}

int calprofile_gtf(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *gtffile, int format, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fgtffile=File_Open(gtffile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char));
    char *strand = malloc(2); int pstrand = 2; //.
    char *geneid = malloc(100*sizeof(char));
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilestep); //*3
    int bodysplitN = ceil(1.0/bodyprofilestep);
    assert(splitN>0 && bodysplitN>0);
    int Total_splitN = 0;
    if(profilemode == 0) {
        Total_splitN = splitN*2 + bodysplitN;
    }else{
        Total_splitN = splitN*2;
    }
    //aver
    int Tsize = 4;
    double *finalstats = calloc(Total_splitN*Tsize, sizeof(double));
    int *finalcounts = calloc(Total_splitN*Tsize, sizeof(int));
    int start=0, end=0, middle = 0, processtart = 0, processend = 0;
    char* storebuffer_c = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_cg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_chg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_chh = malloc(sizeof(char)*MAX_BUFF_PRINT);
    uint32_t Nprocess = 0;
    while(fgets(PerLine,200,Fgtffile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
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
        }
        if(profilemode == 0){ // gene flanks mode
            if(start > upstream) processtart = start-upstream;
            else processtart = 0;
            processend = end+downstream;
            profile_print_array_buffer(finalstats, finalcounts, fp, chrom, processtart, start, end, processend, splitN, bodysplitN, profilemovestep, bodyprofilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
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
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
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
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
        }else if(profilemode == 3){ //center mode
            middle = (int)((start+end)/2);
            if(middle>upstream) processtart = middle - upstream;
            else processtart = 0;
            processend = middle + downstream;
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, splitN*2, profilemovestep, "weighted", pstrand, geneid, context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
        }
    }

    //print last buffer gtf
    if(stored_buffer>0){
        fprintf(stderr, "print last buffer matrix\n");
        fprintf(outfileF_c,"%s",storebuffer_c);
        fprintf(outfileF_cg,"%s",storebuffer_cg);
        fprintf(outfileF_chg,"%s",storebuffer_chg);
        fprintf(outfileF_chh,"%s",storebuffer_chh);
        free(storebuffer_c); free(storebuffer_cg); free(storebuffer_chg); free(storebuffer_chh);
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
    bwClose(fp);
    bwCleanup();
}

int calprofile(char *inbmfile, int upstream, int downstream, double profilestep, double profilemovestep, double bodyprofilestep, double bodyprofilemovestep, char *bedfile, uint8_t context, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, FILE* outfileF_aver, int profilemode, int matrixX){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char));
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilestep); //*3
    int bodysplitN = ceil(1.0/bodyprofilestep);
    assert(splitN>0 && bodysplitN>0);
    int Total_splitN = 0;
    if(profilemode == 0) {
        Total_splitN = splitN*2 + bodysplitN;
    }else{
        Total_splitN = splitN*2;
    }
    //aver
    int Tsize = 4;
    double *finalstats = calloc(Total_splitN*Tsize, sizeof(double));
    int *finalcounts = calloc(Total_splitN*Tsize, sizeof(int));
    int start=0, end=0, middle = 0, processtart = 0, processend = 0;
    char* storebuffer_c = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_cg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_chg = malloc(sizeof(char)*MAX_BUFF_PRINT);
    char* storebuffer_chh = malloc(sizeof(char)*MAX_BUFF_PRINT);
    uint16_t Nprocess = 0;
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
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
        if(profilemode == 0){ // gene flanks mode
            if(start > upstream) processtart = start-upstream;
            else processtart = 0;
            processend = end+downstream;
            profile_print_array_buffer(finalstats, finalcounts, fp, chrom, processtart, start, end, processend, splitN, bodysplitN, profilemovestep, bodyprofilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
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
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
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
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
        }else if(profilemode == 3){ //center mode
            middle = (int)((start+end)/2);
            if(middle>upstream) processtart = middle - upstream;
            else processtart = 0;
            processend = middle + downstream;
            profile_print_array_buffer1(finalstats, finalcounts, fp, chrom, processtart, processend, Total_splitN, profilemovestep, "weighted", pstrand, "", context, storebuffer_c, storebuffer_cg, storebuffer_chg, storebuffer_chh, matrixX);
            if(stored_buffer>MAX_BUFF_PRINT/10000){
                fprintf(outfileF_c,"%s",storebuffer_c);
                storebuffer_c[0] = '\0';
                fprintf(outfileF_cg,"%s",storebuffer_cg);
                storebuffer_cg[0] = '\0';
                fprintf(outfileF_chg,"%s",storebuffer_chg);
                storebuffer_chg[0] = '\0';
                fprintf(outfileF_chh,"%s",storebuffer_chh);
                storebuffer_chh[0] = '\0';
                stored_buffer = 0;
            }
        }
    }

    //print last buffer
    if(stored_buffer>0){
        fprintf(stderr, "print last buffer matrix\n");
        fprintf(outfileF_c,"%s\n",storebuffer_c);
        fprintf(outfileF_cg,"%s",storebuffer_cg);
        fprintf(outfileF_chg,"%s",storebuffer_chg);
        fprintf(outfileF_chh,"%s",storebuffer_chh);
        free(storebuffer_c); free(storebuffer_cg); free(storebuffer_chg); free(storebuffer_chh);
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
    free(print_context);
    free(finalstats); free(finalcounts);

    bwClose(fp);
    bwCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
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

void calbody_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF){
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, end-start, method, pstrand);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    getcontext(context, print_context);
    if(stats) {
        if(format == 1 || format == 2){
            if(context >= 4){
                if(!isnan(stats[0])) {
                    fprintf(outfileF, "C\t%s\t%f\t%s\n", bodycase, stats[0], geneid);
                }
                if(!isnan(stats[1])) {
                    fprintf(outfileF, "CG\t%s\t%f\t%s\n", bodycase, stats[1], geneid);
                }
                if(!isnan(stats[2])) {
                    fprintf(outfileF, "CHG\t%s\t%f\t%s\n", bodycase, stats[2], geneid);
                }
                if(!isnan(stats[3])) {
                    fprintf(outfileF, "CHH\t%s\t%f\t%s\n", bodycase, stats[3], geneid);
                }
            }else if(!isnan(stats[context])) {
                fprintf(outfileF, "%s\t%s\t%f\t%s\n", print_context, bodycase, stats[context], geneid);
            }
        }else{
            if(context >= 4) {
                if(!isnan(stats[0])) {
                    fprintf(outfileF, "C\t%s\t%f\t%s:%d-%d\n", bodycase, stats[0], chrom, start, end);
                }
                if(!isnan(stats[1])) {
                    fprintf(outfileF, "CG\t%s\t%f\t%s:%d-%d\n", bodycase, stats[1], chrom, start, end);
                }
                if(!isnan(stats[2])) {
                    fprintf(outfileF, "CHG\t%s\t%f\t%s:%d-%d\n", bodycase, stats[2], chrom, start, end);
                }
                if(!isnan(stats[3])) {
                    fprintf(outfileF, "CHH\t%s\t%f\t%s:%d-%d\n", bodycase, stats[3], chrom, start, end);
                }
            }else if(!isnan(stats[context])) {
                fprintf(outfileF, "%s\t%s\t%f\t%s:%d-%d\n", print_context, bodycase, stats[context], chrom, start, end);
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

void calregion_weighted_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* bodycase, char* strand, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh, uint16_t *countC, uint16_t *countCT){
    Sregionstats_array_count(fp, chrom, start, end, splitN, end-start, method, pstrand, countC, countCT);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    char* print_geneid = malloc(sizeof(char)*20);
    if(format == 1 || format == 2) {
        strcpy(print_geneid, "\t");
        strcat(print_geneid, geneid);
    }else{
        strcpy(print_geneid, "");
    }
    getcontext(context, print_context);
    if(countCT && countC) {
        if(format == 1 || format == 2 || format == 0){ //bed gtf etc
            if(context >= 4){
                if(countCT[0]>0) { // C
                    fprintf(outfileF_c, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[0], countCT[0], print_geneid);
                }
                if(countCT[1]>0) { // CG
                    fprintf(outfileF_cg, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[1], countCT[1], print_geneid);
                }
                if(countCT[2]>0) {
                    fprintf(outfileF_chg, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[2], countCT[2], print_geneid);
                }
                if(countCT[3]>0) {
                    fprintf(outfileF_chh, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[3], countCT[3], print_geneid);
                }
            }else if(countCT[context]>0) { // print_context
                if(context == 0) fprintf(outfileF_c, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[context], countCT[context], print_geneid);
                else if(context == 1) fprintf(outfileF_cg, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[context], countCT[context], print_geneid);
                else if(context == 2) fprintf(outfileF_chg, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[context], countCT[context], print_geneid);
                else if(context == 3) fprintf(outfileF_chh, "%s\t%d\t%c\t%d\t%d%s\n", chrom, start, strand[0], countC[context], countCT[context], print_geneid);
            }
        }

    }
    free(print_context);
    free(print_geneid);
}

void calregion_print(bigWigFile_t *fp, char* chrom, int start, int end, int splitN, char* method, int pstrand, int format, char* geneid, uint8_t context, char* strand, FILE* outfileF){
    double *stats = Sregionstats_array(fp, chrom, start, end, splitN, end-start, method, pstrand);
    //int i = 0;
    char* print_context = malloc(sizeof(char)*5);
    getcontext(context, print_context);
    if(stats) {
        if(format == 1 || format == 2){ // gtf gff
            if(context >= 4){
                if(!isnan(stats[0])) { //C
                    fprintf(outfileF, "%s\t%d\t%c\t%f\t%s\n", chrom, start, strand[0], stats[0], geneid);
                }
                if(!isnan(stats[1])) { //CG
                    fprintf(outfileF, "%s\t%d\t%c\t%f\t%s\n", chrom, start, strand[0], stats[1], geneid);
                }
                if(!isnan(stats[2])) { //CHG
                    fprintf(outfileF, "%s\t%d\t%c\t%f\t%s\n", chrom, start, strand[0], stats[2], geneid);
                }
                if(!isnan(stats[3])) { //CHH
                    fprintf(outfileF, "%s\t%d\t%c\t%f\t%s\n", chrom, start, strand[0], stats[3], geneid);
                }
            }else if(!isnan(stats[context])) { //print_context
                fprintf(outfileF, "%s\t%d\t%c\t%f\t%s\n", chrom, start, strand[0], stats[context], geneid);
            }
        }else{ //bed
            if(context >= 4){
                if(!isnan(stats[0])) { //C
                    fprintf(outfileF, "%s\t%d\t%c\t%f\n", chrom, start, strand[0], stats[0]);
                }
                if(!isnan(stats[1])) { //CG
                    fprintf(outfileF, "%s\t%d\t%c\t%f\n", chrom, start, strand[0], stats[1]);
                }
                if(!isnan(stats[2])) { //CHG
                    fprintf(outfileF, "%s\t%d\t%c\t%f\n", chrom, start, strand[0], stats[2]);
                }
                if(!isnan(stats[3])) { //CHH
                    fprintf(outfileF, "%s\t%d\t%c\t%f\n", chrom, start, strand[0], stats[3]);
                }
            }else if(!isnan(stats[context])) { //print_context
                fprintf(outfileF, "%s\t%d\t%c\t%f\n", chrom, start, strand[0], stats[context]);
            }
        }
    }
    free(print_context);
}

int calregionstats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, Tsize = 4; //, i =0
    char *geneid = malloc(100*sizeof(char));
    uint16_t *countC = malloc(sizeof(uint16_t)*splitN*Tsize);
    uint16_t *countCT = malloc(sizeof(uint16_t)*splitN*Tsize);
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 0){
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        }else if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
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

    bwClose(fp);
    bwCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
    free(countC); free(countCT);
    return 0;
}


int calregionstats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
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
                calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, 0, "", context, "", strandstr, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
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

    bwClose(fp);
    bwCleanup();
    //free(chrom);
    free(countC); free(countCT);
    return 0;
}

int calbodystats_file(char *inbmfile, char *method, char *bedfile, int format, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, Tsize = 4; //, i =0
    char *geneid = malloc(100*sizeof(char));
    int upstream = 0, downstream = 0;
    uint16_t *countC = malloc(sizeof(uint16_t)*splitN*Tsize);
    uint16_t *countCT = malloc(sizeof(uint16_t)*splitN*Tsize);
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 0){
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        }else if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
            delete_char2(geneid, '"', ';');
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
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
        upstream = start-regionextend>=0?start - regionextend:0;
        downstream = end + regionextend;
        calbody_print(fp, chrom, upstream, start, splitN, method, pstrand, format, geneid, context, "up", strand, outfileF);
        calbody_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, "body", strand, outfileF);
        calbody_print(fp, chrom, end, downstream, splitN, method, pstrand, format, geneid, context, "down", strand, outfileF);
        if(printcoverage > 0){
            calregion_weighted_print(fp, chrom, start, end, splitN, method, pstrand, format, geneid, context, "body", strand, outfileF_c, outfileF_cg, outfileF_chg, outfileF_chh, countC, countCT);
        }
    }

    bwClose(fp);
    bwCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
    free(countC); free(countCT);
    return 0;
}


int calbodystats(char *inbmfile, char *method, char *region, uint8_t pstrand, uint8_t context, FILE* outfileF, FILE* outfileF_c, FILE* outfileF_cg, FILE* outfileF_chg, FILE* outfileF_chh){
    //open mbw file
    uint32_t type = BMtype(inbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(inbmfile, NULL, "r");
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
    uint16_t *countC = malloc(sizeof(uint16_t)*splitN*Tsize);
    uint16_t *countCT = malloc(sizeof(uint16_t)*splitN*Tsize);
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

    bwClose(fp);
    bwCleanup();
    //free(chrom);
    free(countC); free(countCT);
    return 0;
}

int main_view_all(bigWigFile_t *ifp, FILE* outfileF, char *outformat, bigWigFile_t *ofp, char* filterchrom){

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
        //fprintf(stderr, "process %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            sprintf(region, "%s:%d-%d", chrom, start, end);
            //fprintf(stderr, "ccx %s\n", region);
            printL = main_view_mbw(ifp, region, outfileF, outformat, ofp, chromsUse, starts, ends, values, coverages, strands, contexts, 
                entryid);
            if(strcmp(outformat, "mbw") == 0) {
                if(start == 0 && printL>0){
                    int response = bwAddIntervals(ofp, chromsUse, starts, ends, values, coverages, strands, contexts, 
                    entryid, printL);
                    if(response) {
                        fprintf(stderr, "mbwAddIntervals 0\n");
                        return -1;
                    }
                    printL = 0;
                }else if(printL>0){
                    if(bwAppendIntervals(ofp, starts, ends, values, coverages, strands, contexts, entryid, printL)) {
                        fprintf(stderr, "mbwAppendIntervals 2\n");
                        return -1;
                    }
                    printL = 0;
                }
            }
            start += SEGlen;
            end += SEGlen;
        }
    }
    free(region);
    for(i =0; i < MAX_LINE_PRINT; i++){
        free(chromsUse[i]); free(entryid[i]);
    }
    free(chromsUse); free(entryid); free(starts);
    free(ends); free(values); free(coverages); free(strands); free(contexts);
}

int main_view_mbw(bigWigFile_t *ifp, char *region, FILE* outfileF, char *outformat, bigWigFile_t *ofp, \
    char** chromsUse, uint32_t* starts, uint32_t* ends, float* values, uint16_t* coverages, uint8_t* strands, \
    uint8_t* contexts, char** entryid){
    // read. test/example_output.bw
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;
    bwOverlappingIntervals_t *o;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", region, ifp->type);
    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *tempstore = malloc(sizeof(char)*10000000);
    char *tempchar = malloc(20);
    int Nprint = 0;
    //char *strand = malloc(100*sizeof(char)); int strand;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        //strand = strtok(NULL,",:-");
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        o = bwGetOverlappingIntervals(ifp, chrom, start, end+1);
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

                if(strcmp(outformat, "mbw") == 0){
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
                        entryid[Nprint] = o->entryid[j];
                    }
                    Nprint++;
                }else{
                    //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                    //o->strand[i], o->context[i]);   %"PRIu32"
                    sprintf(tempchar, "%s\t%ld", chrom, o->start[j]);
                    strcat(tempstore, tempchar);
                    if(ifp->hdr->version & BM_END) { //ifp->type
                        sprintf(tempchar, "\t%"PRIu32"", o->end[j]);
                        strcat(tempstore, tempchar);
                    }
                    sprintf(tempchar, "\t%f", o->value[j]);
                    strcat(tempstore, tempchar);
                    if(ifp->hdr->version & BM_COVER) {
                        sprintf(tempchar, "\t%"PRIu16"", o->coverage[j]);
                        strcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_STRAND){
                        sprintf(tempchar, "\t%s", strand_str[o->strand[j]]);
                        strcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_CONTEXT) {
                        sprintf(tempchar, "\t%s", context_str[o->context[j]]);
                        strcat(tempstore, tempchar);
                    }
                    if(ifp->hdr->version & BM_ID) {
                        sprintf(tempchar, "\t%s", o->entryid[j]);
                        strcat(tempstore, tempchar);
                    }
                    sprintf(tempchar, "\n");
                    strcat(tempstore, tempchar);
                    Nprint++;
                    if(Nprint>10000) {
                        if(strcmp(outformat, "txt") == 0) fprintf(outfileF,"%s",tempstore);
                        else if(strcmp(outformat, "mbw") == 0) {
                            //mbw out
                        }
                        Nprint = 0;
                    }
                }
            }
        }
    }

    if(Nprint>0) {
        if(strcmp(outformat, "txt") == 0) {
            fprintf(outfileF,"%s",tempstore);
        } else if(strcmp(outformat, "mbw") == 0) {
            ;//mbw out
        }
    }
    //free(chrom);
    free(tempchar);
    free(tempstore);
    bwDestroyOverlappingIntervals(o);
    return Nprint;

error:
    fprintf(stderr, "No results found!\n");
    return 0;
}

int main_view_file(bigWigFile_t *ifp, char *bedfile, FILE* outfileF, char *outformat, bigWigFile_t *ofp){
    // read. test/example_output.bw
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;
    bwOverlappingIntervals_t *o;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", bedfile, ifp->type);

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    unsigned int j = 0;
    char *tempstore = malloc(sizeof(char)*10000000);
    char *tempchar = malloc(20);
    int Nprint = 0;
    char* region = malloc(sizeof(char)*1000);
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
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

    bwDestroyOverlappingIntervals(o);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    return 1;
}

int main_view(bigWigFile_t *ifp, char *region, FILE* outfileF, char *outformat, bigWigFile_t *ofp){
    // read. test/example_output.bw
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;
    bwOverlappingIntervals_t *o;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", region, ifp->type);
    char *substr= strtok(region, ";");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,";");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *tempstore = malloc(sizeof(char)*10000000);
    char *tempchar = malloc(20);
    int Nprint = 0;
    //char *strand = malloc(100*sizeof(char)); int strand;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ",:-");
        start = atoi(strtok(NULL,",:-"));
        end = atoi(strtok(NULL,",:-")); // + 1;
        //strand = strtok(NULL,",:-");
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        o = bwGetOverlappingIntervals(ifp, chrom, start, end+1);
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
                //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                //o->strand[i], o->context[i]);   %"PRIu32"
                sprintf(tempchar, "%s\t%ld", chrom, o->start[j]);
                strcat(tempstore, tempchar);
                if(ifp->hdr->version & BM_END) { //ifp->type
                    sprintf(tempchar, "\t%"PRIu32"", o->end[j]);
                    strcat(tempstore, tempchar);
                }
                sprintf(tempchar, "\t%f", o->value[j]);
                strcat(tempstore, tempchar);
                if(ifp->hdr->version & BM_COVER) {
                    sprintf(tempchar, "\t%"PRIu16"", o->coverage[j]);
                    strcat(tempstore, tempchar);
                }
                if(ifp->hdr->version & BM_STRAND){
                    sprintf(tempchar, "\t%s", strand_str[o->strand[j]]);
                    strcat(tempstore, tempchar);
                }
                if(ifp->hdr->version & BM_CONTEXT) {
                    sprintf(tempchar, "\t%s", context_str[o->context[j]]);
                    strcat(tempstore, tempchar);
                }
                if(ifp->hdr->version & BM_ID) {
                    sprintf(tempchar, "\t%s", o->entryid[j]);
                    strcat(tempstore, tempchar);
                }
                sprintf(tempchar, "\n");
                strcat(tempstore, tempchar);
                Nprint++;
                if(Nprint>10000) {
                    if(strcmp(outformat, "txt") == 0) fprintf(outfileF,"%s",tempstore);
                    else if(strcmp(outformat, "mbw") == 0) {
                        //mbw out
                    }
                    Nprint = 0;
                }
            }
        }
    }

    if(Nprint>0) {
        if(strcmp(outformat, "txt") == 0) fprintf(outfileF,"%s",tempstore);
        else if(strcmp(outformat, "mbw") == 0) {
            //mbw out
        }
        Nprint = 0;
    }
    //free(chrom);
    free(tempchar);
    free(tempstore);
    bwDestroyOverlappingIntervals(o);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    return 1;
}

int main_view_bedfile(char *inbmF, char *bedfile, int type, FILE* outfileF, char *outformat, bigWigFile_t *ofp){
    // read. test/example_output.bw
    bigWigFile_t *ifp = NULL;
    ifp = bwOpen(inbmF, NULL, "r");
    ifp->type = type;
    bwOverlappingIntervals_t *o;

    FILE *BedF = File_Open(bedfile, "r");

    if(DEBUG>1) fprintf(stderr, "xxx1211 %s\n", bedfile);
    int i =0, j = 0;
    
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char region[100];
    while (fgets(region,100,BedF)!=0){
        chrom = strtok(region, ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-"));
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d chrom %s %d %d", i, chrom, start, end);
        o = bwGetOverlappingIntervals(ifp, chrom, start, end+1);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                //o->strand[i], o->context[i]);
                fprintf(stderr, "chr1\t%ld\t%ld\t%f\t%ld\t%s\t%s\t%s\n", o->start[j], o->end[j], o->value[j], o->coverage[j],
                strand_str[o->strand[j]], context_str[o->context[j]], o->entryid[j]);
            }
        }
    }

    //free(chrom);
    fclose(BedF);
    bwDestroyOverlappingIntervals(o);
    bwClose(ifp);
    return 0;

error:
    fprintf(stderr, "No results found!\n");
    bwClose(ifp);
    bwCleanup();
    return 1;
}

int bw_overlap_all_mul(char *inbmFs, uint8_t pstrand){
    char *substr= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        substr = strtok(NULL,",");
    }

    bigWigFile_t **ifps = malloc(sizeof(bigWigFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        bigWigFile_t *ifp1 = NULL;
        ifp1 = bwOpen(infiles[i], NULL, "r");
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
            bw_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
            start += SEGlen;
            end += SEGlen;
        }
    }

    
    for(i=0;i<sizeifp;i++){
        bwClose(ifps[i]);
    }
    return 0;
}

int bw_overlap_file_mul(char *inbmFs, char *bedfile){
    /*
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    bigWigFile_t *ifp1 = NULL;
    ifp1 = bwOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    bigWigFile_t *ifp2 = NULL;
    ifp2 = bwOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;
    */


    char *substr= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        substr = strtok(NULL,",");
    }

    bigWigFile_t **ifps = malloc(sizeof(bigWigFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        bigWigFile_t *ifp1 = NULL;
        ifp1 = bwOpen(infiles[i], NULL, "r");
        ifp1->type = ifp1->hdr->version;
        ifps[i] = ifp1;
    }
    
    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int SEGlen = 1000000, nK = 0, j = 0;
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
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
                bw_overlap_mul(ifps, sizeifp, chrom, start+j*SEGlen, start+(j+1)*SEGlen-1, pstrand);
            }
            bw_overlap_mul(ifps, sizeifp, chrom, start+nK*SEGlen, end, pstrand);
        }else
            bw_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
    }

    for(i=0;i<sizeifp;i++){
        bwClose(ifps[i]);
    }
    fclose(Fbedfile);
    free(PerLine); free(chrom); free(strand);
}

int bw_overlap_file(char *inbmF1, char *inbmF2, char *bedfile){
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    bigWigFile_t *ifp1 = NULL;
    ifp1 = bwOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    bigWigFile_t *ifp2 = NULL;
    ifp2 = bwOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    
    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    //int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        bw_overlap(ifp1, ifp2, chrom, start, end, pstrand);
    }

    bwClose(ifp1);
    bwClose(ifp2);
    fclose(Fbedfile);
    free(PerLine); free(chrom); free(strand);
}

int bw_overlap_region_mul(char *inbmFs, char *region, uint8_t pstrand){
    //open file1
    char *subinbm= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (subinbm != NULL) {
        strcpy(infiles[sizeifp++], subinbm);
        subinbm = strtok(NULL,",");
    }

    bigWigFile_t **ifps = malloc(sizeof(bigWigFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        uint32_t type1 = BMtype(infiles[i], NULL);
        bigWigFile_t *ifp1 = NULL;
        ifp1 = bwOpen(infiles[i], NULL, "r");
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
                bw_overlap_mul(ifps, sizeifp, chrom, start+j*SEGlen, start+(j+1)*SEGlen-1, pstrand);
            }
            bw_overlap_mul(ifps, sizeifp, chrom, start+nK*SEGlen, end, pstrand);
        }else
            bw_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand);
    }

    for(i=0;i<sizeifp;i++){
        bwClose(ifps[i]);
    }
    return 0;
}

int bw_overlap_region(char *inbmF1, char *inbmF2, char *region, uint8_t pstrand){
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    bigWigFile_t *ifp1 = NULL;
    ifp1 = bwOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    bigWigFile_t *ifp2 = NULL;
    ifp2 = bwOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    char *substr= strtok(region, ",");
    char regions[1000][200] = {""};
    int slen = 0, i =0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,",");
    }

    char *chrom = malloc(100*sizeof(char));
    int start=0, end=0;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-")); // + 1;
        bw_overlap(ifp1, ifp2, chrom, start, end, pstrand);
    }

    bwClose(ifp1);
    bwClose(ifp2);
    return 0;
}

int bw_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand){
    if(n1<1 || n2<1){
        return -1;
    }
    //open file1
    uint32_t type1 = BMtype(inbmF1, NULL);
    bigWigFile_t *ifp1 = NULL;
    ifp1 = bwOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    uint32_t type2 = BMtype(inbmF2, NULL);
    bigWigFile_t *ifp2 = NULL;
    ifp2 = bwOpen(inbmF2, NULL, "r");
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
            bw_overlap(ifp1, ifp2, chrom, start, end, pstrand);
            start += SEGlen;
            end += SEGlen;
        }
    }

    bwClose(ifp1);
    bwClose(ifp2);
    return 0;
}

int bw_overlap(bigWigFile_t *ifp1, bigWigFile_t *ifp2, char *chrom, int start, int end, uint8_t strand){
    bwOverlappingIntervals_t *o1;
    bwOverlappingIntervals_t *o2;

    int slen = 1, i =0, j = 0, k = 0, lociK = 0;
    int* countM = malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    for(i=0;i<slen; i++){
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        o1 = bwGetOverlappingIntervals(ifp1, chrom, start, end+1);
        o2 = bwGetOverlappingIntervals(ifp2, chrom, start, end+1);
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
    }

    //free(chrom);
    bwDestroyOverlappingIntervals(o1);
    bwDestroyOverlappingIntervals(o2);
    free(countM);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    return 1;
}

int bw_overlap_mul(bigWigFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand){
    bwOverlappingIntervals_t *o1;
    fprintf(stderr, "process region %s %d %d\n", chrom, start, end);
    int slen = 1, i =0, j = 0;
    int* countM = malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    
    //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
    if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
    int total = 0, loci = 0;
    for(i=0;i<sizeifp;i++){
        o1 = bwGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
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
    }
    if(total == 0){
        bwDestroyOverlappingIntervals(o1);
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
        o1 = bwGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
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
    }

    for(i=0;i<end-start+1;i++){
        if(countM[i] == sizeifp){
            printf("%s\n",printmr[i]);
        }
    }

    //free(chrom);
    bwDestroyOverlappingIntervals(o1);
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

void bwPrintHdr(bigWigFile_t *bw) {
    uint64_t i;
    int64_t i64;
    //fprintf(stderr, "Version:    %"PRIu16"\n", bw->hdr->version);
    if(bw->hdr->version & BM_END) fprintf(stderr, "BM_END:    yes\n");
    else fprintf(stderr, "BM_END:    no\n");
    if(bw->hdr->version & BM_COVER) fprintf(stderr, "BM_COVER:    yes\n");
    else fprintf(stderr, "BM_COVER:    no\n");
    if(bw->hdr->version & BM_CONTEXT) fprintf(stderr, "BM_CONTEXT:    yes\n");
    else fprintf(stderr, "BM_CONTEXT:    no\n");
    if(bw->hdr->version & BM_STRAND) fprintf(stderr, "BM_STRAND:    yes\n");
    else fprintf(stderr, "BM_STRAND:    no\n");
    if(bw->hdr->version & BM_ID) fprintf(stderr, "BM_ID:    yes\n");
    else fprintf(stderr, "BM_ID:    no\n");
    fprintf(stderr, "Levels:     %"PRIu16"\n", bw->hdr->nLevels);
    //fprintf(stderr, "ctOffset:   0x%"PRIx64"\n", bw->hdr->ctOffset);
    //fprintf(stderr, "dataOffset: 0x%"PRIx64"\n", bw->hdr->dataOffset);
    fprintf(stderr, "indexOffset:        0x%"PRIx64"\n", bw->hdr->indexOffset);
    //fprintf(stderr, "sqlOffset:  0x%"PRIx64"\n", bw->hdr->sqlOffset);
    //fprintf(stderr, "summaryOffset:      0x%"PRIx64"\n", bw->hdr->summaryOffset);
    fprintf(stderr, "bufSize:    %"PRIu32"\n", bw->hdr->bufSize);
    fprintf(stderr, "extensionOffset:    0x%"PRIx64"\n", bw->hdr->extensionOffset);

    if(bw->hdr->nLevels) {
        fprintf(stderr, "	i	level	data	index\n");
    }
    for(i=0; i<bw->hdr->nLevels; i++) {
        fprintf(stderr, "\t%"PRIu64"\t%"PRIu32"\t%"PRIx64"\t%"PRIx64"\n", i, bw->hdr->zoomHdrs->level[i], bw->hdr->zoomHdrs->dataOffset[i], bw->hdr->zoomHdrs->indexOffset[i]);
    }

    fprintf(stderr, "nBasesCovered:      %"PRIu64"\n", bw->hdr->nBasesCovered);
    fprintf(stderr, "minVal:     %f\n", bw->hdr->minVal);
    fprintf(stderr, "maxVal:     %f\n", bw->hdr->maxVal);
    //fprintf(stderr, "sumData:    %f\n", bw->hdr->sumData);
    //fprintf(stderr, "sumSquared: %f\n", bw->hdr->sumSquared);

    //Chromosome idx/name/length
    if(bw->cl) {
        fprintf(stderr, "Chromosome List\n");
        fprintf(stderr, "  idx\tChrom\tLength (bases)\n");
        for(i64=0; i64<bw->cl->nKeys; i64++) {
            fprintf(stderr, "  %"PRIu64"\t%s\t%"PRIu32"\n", i64, bw->cl->chrom[i64], bw->cl->len[i64]);
        }
    }
}

void bwPrintIndexNode(bwRTreeNode_t *node, int level) {
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
            bwPrintIndexNode(node->x.child[i], level+1);
        }
    }
}

void bwPrintIndexTree(bigWigFile_t *fp) {
    printf("\nIndex tree:\n");
    printf("nItems:\t%"PRIu64"\n", fp->idx->nItems);
    printf("chrIdxStart:\t%"PRIu32"\n", fp->idx->chrIdxStart);
    printf("baseStart:\t%"PRIu32"\n", fp->idx->baseStart);
    printf("chrIdxEnd:\t%"PRIu32"\n", fp->idx->chrIdxEnd);
    printf("baseEnd:\t%"PRIu32"\n", fp->idx->baseEnd);
    printf("idxSize:\t%"PRIu64"\n", fp->idx->idxSize);
    printf("  level\tchrIdxStart\tbaseStart\tchrIdxEnd\tbaseEnd\tchild\tsize\n");
    bwPrintIndexNode(fp->idx->root, 0);
}