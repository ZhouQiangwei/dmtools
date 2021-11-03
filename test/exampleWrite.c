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
* ./exampleWrite mr2bm -c ~/practice/Genome/hg38/hg38.chr.fa.len -E -C -m test.f -S --Cx -o test.bm -r chr1:0-100,chr1:16766-16890
*/
#include "bigWig.h"
#include <string.h>
#include <stdlib.h>

FILE* File_Open(const char* File_Name,const char* Mode);
char *strand_str[] = {"+", "-", "."};
char *context_str[] = {"C", "CG", "CHG", "CHH"};
int main_view(char *inbmF, char *region);
int main_view_bedfile(char *inbmF, char *bedfile, int type);

#define MAX_LINE_PRINT 100000
int main(int argc, char *argv[]) {
    bigWigFile_t *fp = NULL;
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
    uint32_t i;
    uint32_t write_type = 0x8000;
    char methfile[100]; char *outbmfile = malloc(100);
    char mode[10]; char *region = malloc(1000);
    if(argc>2){
        strcpy(mode, argv[1]);
    }else{
        fprintf(stderr, "Please define mode!!!\n");
    }
    for(i=0; i< argc; i++){
        if(strcmp(argv[i], "-c") == 0){
            strcpy(chromlenf, argv[++i]);
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
        }else if(strcmp(argv[i], "-o") == 0){
            strcpy(outbmfile, argv[++i]);
        }else if(strcmp(argv[i], "-i") == 0){
            strcpy(outbmfile, argv[++i]);
        }else if(strcmp(argv[i], "-r") == 0){
            strcpy(region, argv[++i]);
        }
    }
    if(DEBUG>1) printf("CCCC----- %d %d\n", write_type, sizeof(char*));
    FILE* ChromF=File_Open(chromlenf,"r");
    char *PerLine = malloc(200);
    int printL = 0; char *chrom = malloc(50); int chrlen;
    while(fgets(PerLine,200,ChromF)!=NULL){
        sscanf(PerLine, "%s%d", chrom, &chrLens[printL]);
        chroms[printL++] = strdup(chrom);
    }
    fclose(ChromF);

    printf("[Mode] ------------- %s\n", mode);
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
    if(bwCreateHdr(fp, 10)) goto error;

    //Create the chromosome lists
    fp->cl = bwCreateChromList(chroms, chrLens, printL); //2
    if(!fp->cl) goto error;

    //Write the header
    if(bwWriteHdr(fp)) goto error;
    
    //Some example methlevel
    if(DEBUG>1) fprintf(stderr, "====HHH type %d\n", fp->type);
    char *old_chrom = malloc(50); uint16_t coverC;
    char *strand = malloc(2), *context = malloc(10);
    FILE *methF = File_Open(methfile, "r");
    fgets(PerLine,200,methF); // remove header #
    printL = 0;
    while(fgets(PerLine,200,methF)!=0){
        sscanf(PerLine, "%s%d%s%s%d%d%f", chrom, &starts[printL], strand, context, &coverC, &coverages[printL], &values[printL]);
        chromsUse[printL] = strdup(chrom);
        ends[printL] = starts[printL]+1;
        if(strcmp(old_chrom, chrom)!=0){
            printL = 0;
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
            if(DEBUG>1) fprintf(stderr,"XXX %s %d %d %f %d %d %d\n", chromsUse[0], starts[0], ends[0], values[0], coverages[0], strands[0], context[0]);
            int response = bwAddIntervals(fp, chromsUse, starts, ends, values, coverages, strands, contexts, 
            entryid, 1);
            if(response) goto error;
            if(DEBUG>1) fprintf(stderr, "=response= %d %d %d %d %d\n", response, sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(float));
            strcpy(old_chrom, chrom);
            printL++;
        }else{
            sscanf(PerLine, "%s%d%s%s%d%d%f", chrom, &starts[printL], strand, context, &coverC, &coverages[printL], &values[printL]);
            chromsUse[printL] = strdup(chrom);
            ends[printL] = starts[printL]+1;
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
        if(printL>10000){
            //We can continue appending similarly formatted entries
            //N.B. you can't append a different chromosome (those always go into different
            if(bwAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) goto error;
            printL = 0;
        }
    } // end read me file
    if(printL > 0) {
        if(bwAppendIntervals(fp, starts+1, ends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) goto error;
        printL = 0; 
    }
    
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
    if(DEBUG>0) fprintf(stderr, "bm view1111 ----- \n");
    bwClose(fp);
    if(DEBUG>0) fprintf(stderr, "bm view22222 ---===--- \n");
    bwCleanup();
    /*
    * free memory from malloc
    */
   
    for(i =0; i < MAX_LINE_PRINT; i++){
        free(chroms[i]); free(chromsUse[i]); free(entryid[i]);
    }
    free(chroms); free(chromsUse); free(entryid); 

    free(chrLens); free(starts);
    free(ends); free(values); free(coverages); free(strands); free(contexts);
    free(chrom); free(old_chrom);
    free(strand); free(context); free(PerLine);

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
    if(DEBUG>0) fprintf(stderr, "bm view\n");
    //char region[] = "chr1:0-100,chr1:16766-16830";
    main_view(outbmfile, region);
    free(region);

    //char bedfile[] = "bedfile";
    //main_view_bedfile("test/example_output.bw", bedfile, 4);

    free(chromlenf); free(outbmfile); 
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    bwClose(fp);
    bwCleanup();
    return 1;
}

int main_view(char *inbmF, char *region){
    // read. test/example_output.bw
    uint32_t type = BMtype(inbmF, NULL);
    bigWigFile_t *ifp = NULL;
    ifp = bwOpen(inbmF, NULL, "r");
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d\n", ifp->type);
    ifp->type = type;
    bwOverlappingIntervals_t *o;

    if(DEBUG>1) fprintf(stderr, "xxx111-------- %s %d\n", region, ifp->type);
    char *substr= strtok(region, ",");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,",");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-"));
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        o = bwGetOverlappingIntervals(ifp, chrom, start, end);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                //o->strand[i], o->context[i]);   %"PRIu32"
                printf("--==--==-- %s\t%ld", chrom, o->start[j]);
                if(ifp->type & BM_END)
                    printf("\t%"PRIu32"", o->end[j]);
                printf("\t%f", o->value[j]);
                if(ifp->type & BM_COVER)
                    printf("\t%"PRIu16"", o->coverage[j]);
                if(ifp->type & BM_STRAND)
                    printf("\t%s", strand_str[o->strand[j]]);
                if(ifp->type & BM_CONTEXT)
                    printf("\t%s", context_str[o->context[j]]);
                if(ifp->type & BM_ID)
                    printf("\t%s", o->entryid[j]);
                printf("\n");
            }
        }
    }

    //free(chrom);
    bwDestroyOverlappingIntervals(o);
    bwClose(ifp);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    bwClose(ifp);
    bwCleanup();
    return 1;
}

int main_view_bedfile(char *inbmF, char *bedfile, int type){
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
        o = bwGetOverlappingIntervals(ifp, chrom, start, end);
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
    fprintf(stderr, "Received an error somewhere!\n");
    bwClose(ifp);
    bwCleanup();
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
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}