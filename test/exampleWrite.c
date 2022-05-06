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

#define MAX_LINE_PRINT 1000000
const char* Help_String="Command Format :  bmtools <mode> [opnions]\n"
		"\nUsage:\n"
        "\tmode                  mr2mbw/ view/ overlap/ regionstats/ chromstats.\n"
        "\t [mr2mbw] mode paramaters, required\n"
		"\t-g                    chromosome size file.\n"
        "\t-m                    methratio file\n"
        "\t-o                    output mbigwig file\n"
        "\t [mr2mbw] mode paramaters, options\n"
        "\t-C                    coverage\n"
        "\t-S                    strand\n"
        "\t--Cx                  context\n"
        "\t--Id                  ID\n"
        "\t-E                    end\n"
        "\t--sort Y/N            make chromsize file and meth file in same coordinate, default Y\n"
        "\t--zl                  The maximum number of zoom levels. [1-10]\n"
        "\t-f                    file format. methratio, bedmethyl or bedsimple [default methratio]\n"
        "\t  methratio           chrom start strand context meth_reads cover\n"
        "\t  bedmethyl           chrom start end name * strand * * * coverage meth_reads\n"
        "\t  bedsimple           chrom start end id strand context meth_reads coverage\n"
        "\t--pcontext            CG/CHG/CHH/C, needed when bedmethyl format, default C\n"
        "\t--context             CG/CHG/CHH/ALL, only convert provide context in methratio file or bedsimple, default CG\n"
        "\tNote. meth ratio file must be sorted by chrom and coordinate. ex. sort -k1,1 -k2,2n\n"
        //chrom start end name score strand thickStart thickEnd itemRgb reads_coverage meth_coverage
        "\t [view], [overlap], [regionstats], [profile] mode paramaters, required\n"
        "\t-i                    input mbigwig file\n"
        "\t [view], [overlap], [regionstats] mode paramaters, options\n"
        "\t-r                    region for view, can be seperated by comma.\n"
        "\t--bed                 bed file for view, format: chrom start end [strand].\n"
        "\t [regionstats], [profile] mode paramaters, options\n"
        "\t--gtf                 gtf file for view, format: chrom * * start end * strand * xx geneid.\n"
        "\t--gff                 gff file for view, format: chrom * * start end * strand * xx=geneid.\n"
        "\t [overlap] mode paramaters, required\n"
        "\t-i2                   input mbigwig file2\n"
        "\t [overlap] mode paramaters, required\n"
        "\t--bmfiles             input mbigwig files, seperated by comma.\n"
        "\t [regionstats] [chromstats] mode paramaters, required\n"
        "\t--method              weighted/ mean/ min/ max/ cover/ dev\n"
        "\t [chromstats] mode paramaters, required\n"
        "\t--chromstep           [int] step mean bin size for chromosome region, default: 100000\n"
        "\t--stepmove            [int] step move, default: 50000, if no overlap, please define same as --chromstep\n"

        "\t--pstrand             [0/1/2] strand for calculation, 0 represent '+' positive strand, 1 '-' negative strand, 2 '.' without strand information.\n"
		"\t-h|--help";
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
    char *bmfile2 = malloc(100);
    char mode[10]; char *region = malloc(1000);
    char *inbmfiles = malloc(300);
    int inbm_mul = 0;
    char *method = malloc(100);
    strcpy(method, "weighted");
    int chromstep = 100000;
    int stepoverlap = 50000;
    uint8_t pstrand = 2;
    char *bedfile = malloc(100);
    char *gtffile = malloc(100);
    char *gfffile = malloc(100);

    char *mrformat = malloc(100);
    strcpy(mrformat, "methratio");
    char *pcontext = malloc(10);
    strcpy(pcontext, "C");
    char *filtercontext = malloc(10);
    strcpy(filtercontext, "CG");
    // for gene file
    int upstream = 2000, downstream = 2000;
    double profilestep = 0.02, profilemovestep = 0.01;
    int mcover_cutoff = 4; unsigned long TotalC = 0;
    int zoomlevel = 5;
    char *sortY = malloc(10); strcpy(sortY, "Y");
    if(argc>2){
        strcpy(mode, argv[1]);
        // mr2bam view overlap region
    }else{
        fprintf(stderr, "Please define mode!!!\n");
        fprintf(stderr, "%s\n", Help_String);
        exit(0);
    }
    for(i=0; i< argc; i++){
        if(strcmp(argv[i], "-g") == 0){
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
        }else if(strcmp(argv[i], "--bmfiles") == 0){
            strcpy(inbmfiles, argv[++i]);
            inbm_mul = 1;
        }else if(strcmp(argv[i], "-i2") == 0){
            strcpy(bmfile2, argv[++i]);
        }else if(strcmp(argv[i], "-r") == 0){
            strcpy(region, argv[++i]);
        }else if(strcmp(argv[i], "--method") == 0){
            strcpy(method, argv[++i]);
        }else if(strcmp(argv[i], "--chromstep") == 0){
            chromstep = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--stepmove") == 0){
            stepoverlap = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--pstrand") == 0){
            pstrand = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--bed") == 0){
            strcpy(bedfile, argv[++i]);
        }else if(strcmp(argv[i], "--gtf") == 0){
            strcpy(gtffile, argv[++i]);
        }else if(strcmp(argv[i], "--gff") == 0){
            strcpy(gfffile, argv[++i]);
        }else if(strcmp(argv[i], "-f") == 0){
            strcpy(mrformat, argv[++i]);
        }else if(strcmp(argv[i], "--pcontext") == 0){
            strcpy(pcontext, argv[++i]);
        }else if(strcmp(argv[i], "--sort") == 0){
            strcpy(sortY, argv[++i]);
        }else if(strcmp(argv[i], "--zl") == 0){
            zoomlevel = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--context") == 0){
            strcpy(filtercontext, argv[++i]);
        }
    }

    if(strcmp(mode, "mr2mbw") == 0){
        fprintf(stderr, "mr file format %s\n", mrformat);
        FILE *methF = File_Open(methfile, "r"); 
        char *chrom = malloc(50); char *old_chrom = malloc(50);
        int MAX_CHROM = 10000;
        char *PerLine = malloc(200); char **chromsArray = malloc(sizeof(char*)*MAX_CHROM);
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
        unsigned start=0, end = 0; unsigned coverC,coverage=0; float value=0;
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
                sscanf(PerLine, "%s%u%u%s%*s%s%*s%*s%*s%u%u", chrom, &start, &end, nameid, strand, &coverage, &value);
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
            if(coverage<=mcover_cutoff) continue;
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
                //fprintf(stderr, "\n--222--- %s %d %d\n", chrom, start, coverage);
                chromsUse[printL] = strdup(chrom);
                starts[printL] = start;
                ends[printL] = start+1;
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
        free(chroms); free(chromsUse); free(entryid); 

        free(chrLens); free(starts);
        free(ends); free(values); free(coverages); free(strands); free(contexts);
        free(chrom); free(old_chrom);
        free(strand); free(context); free(PerLine);free(chromlenf); 
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
        printf("[Mode] ------------- View %s\n", region);
        //char region[] = "chr1:0-100,chr1:16766-16830";
        uint32_t type = BMtype(outbmfile, NULL);
        bigWigFile_t *ifp = NULL;
        ifp = bwOpen(outbmfile, NULL, "r");
        ifp->type = ifp->hdr->version;
        if(!region[0] && !bedfile[0]){
            main_view_all(ifp);
        }else if(region[0]){
            main_view(ifp, region);
            free(region);
        }else if(bedfile[0]){
            main_view_file(ifp, bedfile);
            free(bedfile);
        }else{
            fprintf(stderr, "\nplease provide -r or --bed!!!\n");
        }
        free(outbmfile);
        bwClose(ifp);
        return 0;
    }

    //overlap
    if(strcmp(mode, "overlap")==0){
        printf("[Mode] ------------- overlap %s %s\n", outbmfile, bmfile2);
        if(inbm_mul == 1){
            if(!region[0] && !bedfile[0]){
                bw_overlap_all_mul(inbmfiles, pstrand);
            }else if(region[0]){
                bw_overlap_region_mul(inbmfiles, region, pstrand);
                free(region);
            }else if(bedfile[0]){
                bw_overlap_file_mul(inbmfiles, bedfile);
                free(bedfile);
            }else{
                fprintf(stderr, "\nplease provide -r or --bed!!!\n");
            }
            free(inbmfiles); 
        }else{
            if(!region[0] && !bedfile[0]){
                bw_overlap_all(outbmfile, bmfile2, 1, 1, pstrand);
            }else if(region[0]){
                bw_overlap_region(outbmfile, bmfile2, region, pstrand);
                free(region);
            }else if(bedfile[0]){
                bw_overlap_file(outbmfile, bmfile2, bedfile);
                free(bedfile);
            }else{
                fprintf(stderr, "\nplease provide -r or --bed!!!\n");
            }
            free(outbmfile); 
            free(bmfile2); 
        }
        return 0;
    }

    // region
    if(strcmp(mode, "regionstats")==0){
        if(region[0]){
            printf("[Mode] ------------- regionstats %s %s %s\n", outbmfile, region, method);
            calregionstats(outbmfile, method, region, pstrand);
            free(region);
        }else if(bedfile[0]){
            printf("[Mode] ------------- regionstats %s %s %s\n", outbmfile, region, bedfile);
            calregionstats_file(outbmfile, method, bedfile, 0);
            free(bedfile);
        }else if(gtffile[0]){
            printf("[Mode] ------------- regionstats %s %s %s\n", outbmfile, region, gtffile);
            calregionstats_file(outbmfile, method, gtffile, 1);
            free(gtffile);
        }else if(gfffile[0]){
            printf("[Mode] ------------- regionstats %s %s %s\n", outbmfile, region, gfffile);
            calregionstats_file(outbmfile, method, gfffile, 2);
            free(gfffile);
        }else{
            fprintf(stderr, "\nplease provide -r or --bed!!!\n");
        }
        free(outbmfile); 
        return 0;
    }

    // chromstats
    if(strcmp(mode, "chromstats")==0){
        printf("[Mode] ------------- chromstats %s %s\n", outbmfile, method);
        calchromstats(outbmfile, method, chromstep, stepoverlap, pstrand);
        free(outbmfile); 
        free(method);
        return 0;
    }

    //gene profile
    if(strcmp(mode, "profile")==0){
        printf("[Mode] ------------- profile %s\n", outbmfile);
        if(bedfile[0]){
            calprofile(outbmfile, upstream, downstream, profilestep, profilemovestep, bedfile);
            free(bedfile);
        }
        else if(gtffile[0]){
            calprofile_gtf(outbmfile, upstream, downstream, profilestep, profilemovestep, gtffile, 1);
            free(gtffile);
        }else if(gfffile[0]){
            calprofile_gtf(outbmfile, upstream, downstream, profilestep, profilemovestep, gfffile, 2);
            free(gfffile);
        }
        return 0;
    }

    return 1;
error:
    fprintf(stderr, "Received an error in process!\n");
    bwClose(fp);
    bwCleanup();
    return -1;
}

double *Sregionstats(bigWigFile_t *fp, char *chrom, int start, int end, int splitN, uint32_t movestep, char *method, uint8_t strand){
    double *stats = NULL;
    assert(splitN>0);
    int i=0;
    if(strcmp(method, "mean")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, mean, strand);
    }else if(strcmp(method, "weighted")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, weighted, strand);
    }else if(strcmp(method, "dev")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, dev, strand);
    }else if(strcmp(method, "min")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, min, strand);
    }else if(strcmp(method, "max")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, max, strand);
    }else if(strcmp(method, "cover")==0){
        stats = bwStats(fp, chrom, start, end, splitN, movestep, cov, strand);
    }
    return stats;
}

int calchromstats(char *outbmfile, char *method, int chromstep, int stepoverlap, uint8_t strand){
    //open mbw file
    uint32_t type = BMtype(outbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(outbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    int i = 0, j = 0, start = 0, end = chromstep;
    char* region = malloc(sizeof(char)*1000);
    int splitN = 1;
    for(i=0;i<fp->cl->nKeys;i++){
        char* chrom = (char*)fp->cl->chrom[i];
        int len = (int)fp->cl->len[i];
        start = 0; end = chromstep;
        //fprintf(stderr, "CCCC %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            double *stats = Sregionstats(fp, chrom, start, end, splitN, end-start, method, strand);
            if(stats && !isnan(stats[0])) {
                printf("%s:%d-%d\t%f", chrom, start, end, stats[0]);
                if(splitN>1){
                    for(j=1;j<splitN;j++){
                        printf("\t%f", stats[j]);
                    }
                }
                printf("\n");
            }
            start += stepoverlap;
            end += stepoverlap;
        }
    }
    free(region);

    bwClose(fp);
    bwCleanup();
    return 0;
}

int calprofile_gtf(char *outbmfile, int upstream, int downstream, double profilestep, double profilemovestep, char *gtffile, int format){
    //open mbw file
    uint32_t type = BMtype(outbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(outbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fgtffile=File_Open(gtffile,"r");
    char *PerLine = malloc(200);
    int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    char *geneid = malloc(100*sizeof(char));
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilestep)*3;
    assert(splitN>0);
    while(fgets(PerLine,200,Fgtffile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
        }else{
            fprintf(stderr, "\nE: unexpected file format!!!\n");
        }
        fprintf(stderr, "%s %d %d %s %s %s\n", chrom, start, end, strand, geneid);
        if(end-start < splitN) continue;

        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        if(start > upstream) start -= upstream;
        else start = 0;
        end += downstream;
        double *stats = Sregionstats(fp, chrom, start, end, splitN, (int)((end-start)*profilemovestep), "weighted", pstrand);
        if(stats) {
            if(pstrand == 1) {// -
                printf("%s:%d-%d-%s\t%f", chrom, start, end, geneid, stats[splitN-1]);
                if(splitN>1){
                    for(i=splitN-2;i>0;i--){
                        printf("\t%f", stats[i]);
                    }
                }
                printf("\n");
            }else{
                printf("%s:%d-%d-%s\t%f", chrom, start, end, geneid, stats[0]);
                if(splitN>1){
                    for(i=1;i<splitN;i++){
                        printf("\t%f", stats[i]);
                    }
                }
                printf("\n");
            }
        }
    }

    bwClose(fp);
    bwCleanup();
    fclose(Fgtffile);
    free(chrom); free(PerLine); free(strand); free(geneid);
}

int calprofile(char *outbmfile, int upstream, int downstream, double profilestep, double profilemovestep, char *bedfile){
    //open mbw file
    uint32_t type = BMtype(outbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(outbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, i =0;
    splitN = ceil(1.0/profilestep)*3;
    assert(splitN>0);
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        if(end-start < splitN) continue;

        if(strand[0] == '+'){
            pstrand = 0;
        }else if(strand[0] == '-'){
            pstrand = 1;
        }else{
            pstrand = 2;
        }
        if(start > upstream) start -= upstream;
        else start = 0;
        end += downstream;
        double *stats = Sregionstats(fp, chrom, start, end, splitN, (int)((end-start)*profilemovestep), "weighted", pstrand);
        if(stats) {
            if(pstrand == 1) {// -
                printf("%s:%d-%d\t%f", chrom, start, end, stats[splitN-1]);
                if(splitN>1){
                    for(i=splitN-2;i>0;i--){
                        printf("\t%f", stats[i]);
                    }
                }
                printf("\n");
            }else{
                printf("%s:%d-%d\t%f", chrom, start, end, stats[0]);
                if(splitN>1){
                    for(i=1;i<splitN;i++){
                        printf("\t%f", stats[i]);
                    }
                }
                printf("\n");
            }
        }
    }

    bwClose(fp);
    bwCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
}

int calregionstats_file(char *outbmfile, char *method, char *bedfile, int format){
    //open mbw file
    uint32_t type = BMtype(outbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(outbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    FILE* Fbedfile=File_Open(bedfile,"r");
    char *PerLine = malloc(200);
    int printL = 0;
    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    char *strand = malloc(2); int pstrand = 2; //.
    int splitN = 1, i =0;
    char *geneid = malloc(100*sizeof(char));
    while(fgets(PerLine,200,Fbedfile)!=NULL){
        if(PerLine[0] == '#') continue;
        if(format == 0){
            sscanf(PerLine, "%s%d%d%s", chrom, &start, &end, strand);
        }else if(format == 1){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*s%s", chrom, &start, &end, strand, geneid);
        }else if(format == 2){
            sscanf(PerLine, "%s\t%*s\t%*s\t%d\t%d\t%*s\t%s\t%*s\t%*[^=]=%[^;\n\t]", chrom, &start, &end, strand, geneid);
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
        double *stats = Sregionstats(fp, chrom, start, end, splitN, end-start, method, pstrand);
        if(stats) {
            if(format == 1 || format == 2){
                printf("%s:%d-%d-%s\t%f", chrom, start, end, geneid, stats[0]);
            }else{
                printf("%s:%d-%d\t%f", chrom, start, end, stats[0]);
            }
            if(splitN>1){
                for(i=1;i<splitN;i++){
                    printf("\t%f", stats[i]);
                }
            }
            printf("\n");
        }
    }

    bwClose(fp);
    bwCleanup();
    fclose(Fbedfile);
    free(chrom); free(PerLine); free(strand);
    return 0;
}


int calregionstats(char *outbmfile, char *method, char *region, uint8_t pstrand){
    //open mbw file
    uint32_t type = BMtype(outbmfile, NULL);
    bigWigFile_t *fp = NULL;
    fp = bwOpen(outbmfile, NULL, "r");
    fp->type = fp->hdr->version;

    char *substr= strtok(region, ",");
    char regions[1000][200] = {""};
    int slen = 0, i =0, j = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,",");
    }

    char *chrom = malloc(100*sizeof(char)); int start=0, end=0;
    int splitN = 1;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-")); // + 1;
        double *stats = Sregionstats(fp, chrom, start, end, splitN, end-start, method, pstrand);
        if(stats) {
            printf("%s:%d-%d\t%f", chrom, start, end, stats[0]);
            if(splitN>1){
                for(j=1;j<splitN;j++){
                    printf("\t%f", stats[j]);
                }
            }
            printf("\n");
        }
    }

    bwClose(fp);
    bwCleanup();
    //free(chrom);
    return 0;
}

int main_view_all(bigWigFile_t *ifp){

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    int i=0;
    char* region = malloc(sizeof(char)*1000);
    for(i=0;i<ifp->cl->nKeys;i++){
        char* chrom = (char*)ifp->cl->chrom[i];
        int len = (int)ifp->cl->len[i];
        start = 0, end = SEGlen-1;
        //fprintf(stderr, "CCCC %s\t%ld\n", chrom, len);
        while(start<len){
            if(end>len){
                end = len;
            }
            sprintf(region, "%s:%d-%d", chrom, start, end);
            //fprintf(stderr, "ccx %s", region);
            main_view(ifp, region);
            start += SEGlen;
            end += SEGlen;
        }
    }
    free(region);
}

int main_view_file(bigWigFile_t *ifp, char *bedfile){
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
    int j = 0;
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
        o = bwGetOverlappingIntervals(ifp, chrom, start, end+1);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        //fprintf(stderr, "--- version --- %d\n", ifp->hdr->version);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                //o->strand[i], o->context[i]);   %"PRIu32"
                if(pstrand!=2){
                    if(ifp->hdr->version & BM_STRAND){
                        if(o->strand[j] != pstrand){
                            continue;
                        }
                    }
                }
                printf("%s\t%ld", chrom, o->start[j]);
                if(ifp->hdr->version & BM_END) //ifp->type
                    printf("\t%"PRIu32"", o->end[j]);
                printf("\t%f", o->value[j]);
                if(ifp->hdr->version & BM_COVER)
                    printf("\t%"PRIu16"", o->coverage[j]);
                if(ifp->hdr->version & BM_STRAND)
                    printf("\t%s", strand_str[o->strand[j]]);
                if(ifp->hdr->version & BM_CONTEXT)
                    printf("\t%s", context_str[o->context[j]]);
                if(ifp->hdr->version & BM_ID)
                    printf("\t%s", o->entryid[j]);
                printf("\n");
            }
        }
    }

    fclose(Fbedfile);
    free(chrom);
    free(strand);
    free(PerLine);

    bwDestroyOverlappingIntervals(o);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
    return 1;
}

int main_view(bigWigFile_t *ifp, char *region){
    // read. test/example_output.bw
    if(DEBUG>1) fprintf(stderr, "\nifp===-=== %d %d\n", ifp->type, ifp->hdr->version);
    //ifp->type = type;
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
        end = atoi(strtok(NULL,":-")); // + 1;
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        o = bwGetOverlappingIntervals(ifp, chrom, start, end+1);
        if(!o) goto error;
        if(DEBUG>1) fprintf(stderr, "\no->l %ld %ld %d\n", o->l, o->m, ifp->type);
        //fprintf(stderr, "--- version --- %d\n", ifp->hdr->version);
        if(o->l) {
            for(j=0; j<o->l; j++) {
                //fprintf(stderr, "1\t%ld\t%ld\t%f\t%ld\t%d\t%d\n", o->start[i], o->end[i], o->value[i], o->coverage[i],
                //o->strand[i], o->context[i]);   %"PRIu32"
                printf("%s\t%ld", chrom, o->start[j]);
                if(ifp->hdr->version & BM_END) //ifp->type
                    printf("\t%"PRIu32"", o->end[j]);
                printf("\t%f", o->value[j]);
                if(ifp->hdr->version & BM_COVER)
                    printf("\t%"PRIu16"", o->coverage[j]);
                if(ifp->hdr->version & BM_STRAND)
                    printf("\t%s", strand_str[o->strand[j]]);
                if(ifp->hdr->version & BM_CONTEXT)
                    printf("\t%s", context_str[o->context[j]]);
                if(ifp->hdr->version & BM_ID)
                    printf("\t%s", o->entryid[j]);
                printf("\n");
            }
        }
    }

    //free(chrom);
    bwDestroyOverlappingIntervals(o);
    return 0;

error:
    fprintf(stderr, "Received an error somewhere!\n");
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
    fprintf(stderr, "Received an error somewhere!\n");
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
    int printL = 0;
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
    int printL = 0;
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
    char *substr= strtok(region, ",");
    char regions[1000][200] = {""};
    int slen = 0;
    
    while (substr != NULL) {
        strcpy(regions[slen++], substr);
        substr = strtok(NULL,",");
    }

    int SEGlen = 1000000;
    char *chrom = malloc(100*sizeof(char));
    int start=0, end=0, nK = 0, j = 0;
    for(i=0;i<slen; i++){
        chrom = strtok(regions[i], ":-");
        start = atoi(strtok(NULL,":-"));
        end = atoi(strtok(NULL,":-")); // + 1;
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
                        sprintf(tempchar,"\t%d", (int)((double)o1->value[j]*o1->coverage[j]));
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
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}
