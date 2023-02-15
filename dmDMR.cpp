/*
* This tools is used for find dmr with bigmeth file.
* 
* make clean && make && gcc libBigWig.a libBigWig.so test/exampleWrite.c -o exampleWrite
* ./exampleWrite mr2bm -g ~/practice/Genome/hg38/hg38.chr.fa.len -E -C -m test.f -S --Cx -o test.bm -r chr1:0-100,chr1:16766-16890
* g++ test/bmDMR.cpp -o bmDMR -I. -L. -lBigWig -Wl,-rpath /public/home/qwzhou/software_devp/batmeth2/src/bmtools/ -lgsl -lgslcblas -lm -lz
*/
#include "binaMeth.h"
#ifdef __cplusplus
extern "C"
{
#endif
    #include "binaMeth.h"
#ifdef __cplusplus
}
#endif
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_cdf.h>
#include <tr1/cmath>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include "regression.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <iterator>
#include "stddef.h"
#include <map>
//#include <omp.h>

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::ios_base;
using std::cout;
using std::ifstream; 
using std::ofstream; 
using std::floor;
bool geneid = false;

struct PvalLocus {
  std::size_t pos;
  double raw_pval;
  double adjust_pval;
  double combined_pval;
  double corrected_pval;
  size_t position, coverage_factor, meth_factor, coverage_rest, meth_rest;
  std::string chrom, context, sign,name;
};

FILE* File_Open(const char* File_Name,const char* Mode);
const char *strand_str[] = {"+", "-", "."};
const char *context_str[] = {"C", "CG", "CHG", "CHH"};
int bm_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand, vector<PvalLocus> &pvals, unsigned long &NcountC, float Pcutoff, float methdiff, char *filtercontext);
int bm_overlap_all_mul(char *inbmFs, uint8_t pstrand);
int bm_overlap_all_mul_sep(char *inbmFs1, char *inbmFs2, uint8_t pstrand, vector<PvalLocus> &pvals, unsigned long &NcountC, float Pcutoff, float methdiff, char* filtercontext);
int bm_overlap_mul(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand, int Nsample1);
int bm_overlap_mul_calp(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand, int Nsample1, Regression &full_regression, Regression &null_regression, vector<PvalLocus> &pvals, unsigned long &NcountC, float Pcutoff, float methdiff, unsigned long chrom_offset, char* filtercontext);
int bm_overlap(binaMethFile_t *ifp1, binaMethFile_t *ifp2, char *chrom, int start, int end, uint8_t strand,vector<PvalLocus> &loci,
    unsigned long &NcountC, float Pcutoff, float methdiff, unsigned long chrom_offset, char* filtercontext);
void adjust(vector<PvalLocus> &loci, unsigned long NcountC);
void Print_dm_result(const vector<PvalLocus> &pval_loci, ostream &output_encoding, double cutoff, double Pcutoff, double methdiff, string dmr_outfile, 
    bool singleAuto, int mindmc, int mindmcdis, int maxdmrlen, bool geneid);
static double fishers_exact(size_t a, size_t b, size_t c, size_t d);
double loglikratio_test(double null_loglik, double full_loglik);
int checkcomma(char* inbmfiles);

#define MAX_LINE_PRINT 1000000
const char* Help_String="dmDMR [options] -p <prefix of result> -1 [Sample1-methdm,..] -2 [sample2-methdm,..] \n"
                "\nExample:\n"
        		"\tdmDMR -1 s1.dm -2 s2.dm -p prefix --mindmc 5 --minstep 200\n"
                "\nUsage:\n"
        		"\t-p            output file prefix\n"
                "\t-1            sample1 methy dm files, sperate by comma.\n"
                "\t-2            sample2 methy dm files, sperate by comma.\n"
                "\t--mindmc      min dmc sites in dmr region. [default : 4]\n"
                "\t--minstep     min step in bp [default : 100]\n"
                "\t--maxdis      max length of dmr [default : 0]\n"
                "\t--pvalue      pvalue cutoff, default: 0.01\n"
                "\t--fdr         adjust pvalue cutoff default : 1.0\n"
                "\t--methdiff    the cutoff of methylation differention. default: 0.25 [CpG]\n"
                "\t--element     caculate gene or TE etc function elements.\n"
                "\t--context     Context for DM. C/CG/CHG/CHH, [C]\n"
                //"\t-L           predefinded regions or loci.\n"
                "\t-h|--help";
int main(int argc, char *argv[]) {
    binaMethFile_t *fp = NULL;

    uint32_t i;
    uint32_t write_type = 0x8000;
    char methfile[100]; char *inbmfile = (char *)malloc(100);
    char *bmfile2 = (char *)malloc(100);
    char *inbmfiles = (char *)malloc(300);
    int inbm_mul = 0;
    char *method = (char *)malloc(100);
    strcpy(method, "weighted");
    uint8_t pstrand = 2;

    char *mrformat = (char *)malloc(100);
    strcpy(mrformat, "methratio");
    char *filtercontext = (char *)malloc(10);
    strcpy(filtercontext, "C");

    string dmc_outfile;
    string dmr_outfile;
    string prefix = "";    
    double Pcutoff=0.01; double fdr_cutoff = 1; 
    double methdiff=0.25;
    double mdcg = 0.2; double mdchg = 0.1; double mdchh = 0.1;
    bool singleAuto = true;
    int mindmc = 4;
    int mindmcdis = 100;
    int maxdmrlen = 1000000000;
    // for gene file
    if(argc<2){
        fprintf(stderr, "Please define mode!!!\n");
        fprintf(stderr, "%s\n", Help_String);
        exit(0);
    }
    for(i=1; i< argc; i++){
        if(strcmp(argv[i], "--minstep") == 0){
            mindmcdis = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--mindmc") == 0){
            mindmc = atoi(argv[++i]);
        }else if(strcmp(argv[i], "--maxdis") == 0){
            maxdmrlen=atoi(argv[++i]);
        }else if(strcmp(argv[i], "-p") == 0){
            prefix = argv[++i];
            dmc_outfile = prefix + ".dmc";
            dmr_outfile = prefix + ".dmr";
        }else if(strcmp(argv[i], "--Id") == 0){
            write_type |= BM_ID;
        }else if(strcmp(argv[i], "-E") == 0){
            write_type |= BM_END;
        }else if(strcmp(argv[i], "-m") == 0){
            strcpy(methfile, argv[++i]);
        }else if(strcmp(argv[i], "-1") == 0){
            strcpy(inbmfile, argv[++i]);
        }else if(strcmp(argv[i], "--bmfiles") == 0){
            strcpy(inbmfiles, argv[++i]);
            inbm_mul = 1;
        }else if(strcmp(argv[i], "-2") == 0){
            strcpy(bmfile2, argv[++i]);
        }else if(strcmp(argv[i], "--pvalue") == 0){
            Pcutoff = atof(argv[++i]);
        }else if(strcmp(argv[i], "--fdr") == 0){
            fdr_cutoff = atof(argv[++i]);
        }else if(strcmp(argv[i], "--methdiff") == 0){
            methdiff = atof(argv[++i]);
        }else if(strcmp(argv[i], "--context") == 0){
            strcpy(filtercontext, argv[++i]);
        }else if(strcmp(argv[i], "--element") == 0){
            geneid=true;
        }else{
            fprintf(stderr, "Error: unknown paramater %s\n", argv[i]);
            exit(0);
        }
    }

    //overlap
    printf("[Mode] ------------- overlap %s %s\n", inbmfile, bmfile2);
    if(inbm_mul == 1){
        bm_overlap_all_mul(inbmfiles, pstrand);
        free(inbmfiles); 
    }else{
        static vector<PvalLocus> pvals;
        unsigned long NcountC = 0;
        fprintf(stderr, "start overlap and calculate p-value for cytosine site with cutoff %f %f %f\n", Pcutoff, fdr_cutoff, methdiff);
        if(checkcomma(inbmfile) || checkcomma(bmfile2)){
            fprintf(stderr, "with replication, run beta-bionormal test\n");
            bm_overlap_all_mul_sep(inbmfile, bmfile2, pstrand, pvals, NcountC, Pcutoff, methdiff, filtercontext);
        }else{
            fprintf(stderr, "without replication, run p-value test\n");
            bm_overlap_all(inbmfile, bmfile2, 1, 1, pstrand, pvals, NcountC, Pcutoff, methdiff, filtercontext);
        }
        if(pvals.size() == 0){
            fprintf(stderr, "The size of meet pval cutoff is 0, so we exit now.\n");
            exit(0);
        }
        fprintf(stderr, "end overlap and start calculate adjusted p-value %ld, %ld\n", pvals.size(), NcountC);
        ofstream OutFileAdjust;
        adjust(pvals, NcountC);
        fprintf(stderr, "done and print result\n");
        OutFileAdjust.open(dmc_outfile.c_str());
        Print_dm_result(pvals, OutFileAdjust, fdr_cutoff, Pcutoff, methdiff, dmr_outfile, singleAuto, mindmc, mindmcdis, maxdmrlen, geneid);

        free(inbmfile); 
        free(bmfile2); 
    }
    return 0;
}

int checkcomma(char* inbmfiles) {
    for(;*(inbmfiles+1);inbmfiles++){
        if(*inbmfiles == ',') return 1;
    }
    return 0;
}

int bm_overlap_all_mul_sep(char *inbmFs1, char *inbmFs2, uint8_t pstrand, vector<PvalLocus> &pvals, unsigned long &NcountC, float Pcutoff, float methdiff, char* filtercontext){
    int Nsample1 = 0;
    char *substr= strtok(inbmFs1, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;

    Regression full_regression;
    full_regression.design.factor_names.push_back("base");
    full_regression.design.factor_names.push_back("case");

    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        full_regression.design.sample_names.push_back(substr);
        vector<double> matrix_row;
        matrix_row.push_back(1);
        matrix_row.push_back(0);
        full_regression.design.matrix.push_back(vector<double>());
        swap(full_regression.design.matrix.back(), matrix_row);
        substr = strtok(NULL,",");
    }

    Nsample1 = sizeifp;
    substr= strtok(inbmFs2, ",");
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        full_regression.design.sample_names.push_back(substr);
        vector<double> matrix_row;
        matrix_row.push_back(1);
        matrix_row.push_back(1);
        full_regression.design.matrix.push_back(vector<double>());
        swap(full_regression.design.matrix.back(), matrix_row);
        substr = strtok(NULL,",");
    }

    string test_factor_name="case";
    vector<string>::const_iterator test_factor_it =
      std::find(full_regression.design.factor_names.begin(),
                full_regression.design.factor_names.end(), test_factor_name);

    if (test_factor_it == full_regression.design.factor_names.end())
        fprintf(stderr, "Error: %s is not a part of the design specification.", test_factor_name.c_str());
    size_t test_factor = test_factor_it -
                                  full_regression.design.factor_names.begin();

    Regression null_regression;
    null_regression.design = full_regression.design;
    remove_factor(null_regression.design, test_factor);


    binaMethFile_t **ifps = (binaMethFile_t **)malloc(sizeof(binaMethFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
        binaMethFile_t *ifp1 = NULL;
        ifp1 = bmOpen(infiles[i], NULL, "r");
        ifp1->type = ifp1->hdr->version;
        ifps[i] = ifp1;
    }

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    unsigned long chrom_offset = 0;
    for(i=0;i<ifps[0]->cl->nKeys;i++){
        char* chrom = (char*)ifps[0]->cl->chrom[i];
        int len = (int)ifps[0]->cl->len[i];
        start = 0, end = SEGlen-1;
        fprintf(stderr, "process chrom, %s\t%ld\t%d\n", chrom, len, sizeifp);
        while(start<len){
            if(end>len){
                end = len;
            }
            fprintf(stderr, "processing %s:%d-%d\n", chrom, start, end);
            bm_overlap_mul_calp(ifps, sizeifp, chrom, start, end, pstrand, Nsample1, full_regression, null_regression, pvals, NcountC, Pcutoff, methdiff, chrom_offset, filtercontext);
            start += SEGlen;
            end += SEGlen;
        }
        chrom_offset+=len;
    }


    for(i=0;i<sizeifp;i++){
        bmClose(ifps[i]);
    }
    return 0;
}

int bm_overlap_all_mul(char *inbmFs, uint8_t pstrand){
    
    char *substr= strtok(inbmFs, ",");
    char infiles[1000][200] = {""};
    int sizeifp = 0, i = 0;
    
    while (substr != NULL) {
        strcpy(infiles[sizeifp++], substr);
        substr = strtok(NULL,",");
    }

    binaMethFile_t **ifps = (binaMethFile_t **)malloc(sizeof(binaMethFile_t)*sizeifp);

    for(i=0;i<sizeifp;i++){
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
        fprintf(stderr, "process %s\t%ld\t%d\n", chrom, len, sizeifp);
        while(start<len){
            if(end>len){
                end = len;
            }
            fprintf(stderr, "region %s %d %d\n", chrom, start, end);
            bm_overlap_mul(ifps, sizeifp, chrom, start, end, pstrand, 0);
            start += SEGlen;
            end += SEGlen;
        }
    }

    
    for(i=0;i<sizeifp;i++){
        bmClose(ifps[i]);
    }
    return 0;
}

int bm_overlap_all(char *inbmF1, char *inbmF2, int n1, int n2, uint8_t pstrand, vector<PvalLocus> &pvals, unsigned long &NcountC,
    float Pcutoff, float methdiff, char* filtercontext){
    
    if(n1<1 || n2<1){
        return -1;
    }
    //open file1
    binaMethFile_t *ifp1 = NULL;
    ifp1 = bmOpen(inbmF1, NULL, "r");
    ifp1->type = ifp1->hdr->version;
    //file2
    binaMethFile_t *ifp2 = NULL;
    ifp2 = bmOpen(inbmF2, NULL, "r");
    ifp2->type = ifp2->hdr->version;

    int SEGlen = 1000000;
    int start = 0, end = SEGlen-1;
    int i=0; unsigned long chrom_offset = 0;
    for(i=0;i<ifp1->cl->nKeys;i++){
        char* chrom = (char*)ifp1->cl->chrom[i];
        int len = (int)ifp1->cl->len[i];
        start = 0, end = SEGlen-1;
        fprintf(stderr, "process %s\t%ld\n", chrom, len);
        while(start<len) {
            end = start + SEGlen-1;
            if(end>len){
                end = len;
            }
            fprintf(stderr, "region %s %d %d\n", chrom, start, end);
            bm_overlap(ifp1, ifp2, chrom, start, end, pstrand, pvals, NcountC, Pcutoff, methdiff, chrom_offset, filtercontext);
            start += SEGlen;
            end += SEGlen;
        }
        chrom_offset+=len;
    }

    bmClose(ifp1);
    bmClose(ifp2);
    return 0;
}

int bm_overlap(binaMethFile_t *ifp1, binaMethFile_t *ifp2, char *chrom, int start, int end, uint8_t strand,vector<PvalLocus> &pvals,
    unsigned long &NcountC, float Pcutoff, float methdiff, unsigned long chrom_offset, char* filtercontext){
    
    //omp_set_num_threads(10);

    int slen = 1, i =0, j = 0, k = 0, lociK = 0;
    int* countM = (int *)malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    int countC1=0, cover1=0, countC2=0, cover2=0;
    for(i=0;i<slen; i++){
        //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
        if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
        bmOverlappingIntervals_t *o1;
        bmOverlappingIntervals_t *o2;
        o1 = bmGetOverlappingIntervals(ifp1, chrom, start, end+1);
        o2 = bmGetOverlappingIntervals(ifp2, chrom, start, end+1);
        if(!o1 || !o2) {
            fprintf(stderr, "Received an error somewhere!\n");
            return -1;
        };
        //fprintf(stderr, "--- version --- %d\n", ifp1->hdr->version);
        if(o1->l && o2->l) {
            //#pragma omp parallel for
            for(j=0; j<o1->l; j++) {
                if(strand!=2){
                    if(ifp1->hdr->version & BM_STRAND){
                        if(strand != o1->strand[j]){
                            continue;
                        }
                    }
                }
                if(strcmp(filtercontext, "C") != 0){
                    if(ifp1->hdr->version & BM_CONTEXT){
                        if(strcmp(filtercontext, context_str[o1->context[j]]) != 0 ){
                            continue;
                        }
                    }
                }
                //if(o1->start[j] == 2144615) {fprintf(stderr, "1111 %s %d\n", chrom, countM[o1->start[j]-start]);}
                countM[o1->start[j]-start]++;
            }
            //#pragma omp parallel for
            for(j=0; j<o2->l; j++) {
                if(strand!=2){
                    if(ifp2->hdr->version & BM_STRAND){
                        if(strand != o2->strand[j]){
                            continue;
                        }
                    }
                }
                if(strcmp(filtercontext, "C") != 0){
                    if(ifp2->hdr->version & BM_CONTEXT){
                        if(strcmp(filtercontext, context_str[o2->context[j]]) != 0 ){
                            continue;
                        }
                    }
                }
                countM[o2->start[j]-start]++;
                //if(o2->start[j] == 2144615) {fprintf(stderr, "2222 %s %d\n", chrom, countM[o2->start[j]-start]);}
            }
            for(j=0; j<o1->l; j++) {
                if(countM[o1->start[j]-start] == 2){
                    for(k=lociK; k<o2->l; k++) {
                        //if(o1->start[j] == 2144615) {fprintf(stderr, "1111 %s\n", chrom);}
                        //printf("%d %d %d %d\n", j, k, o1->start[j], o2->start[k]);
                        if(o2->start[k]>o1->start[j]) break;
                        if(o1->start[j] != o2->start[k]){
                            continue;
                        }

                        //if(o2->start[k] == 2144615) {fprintf(stderr, "2222 %s\n", chrom);}
                        lociK = k;
                        cover1 = o1->coverage[j];
                        cover2 = o2->coverage[k];
                        countC1 = (int)((double)o1->value[j]*cover1 + 0.5);
                        countC2 = (int)((double)o2->value[k]*cover2 + 0.5);
                        double meth_diff=fabs(double(countC1)/cover1 - double(countC2)/cover2);
                        double pval = 1;
                        if(meth_diff>=methdiff)
                            pval = fishers_exact(countC1,cover1,countC2,cover2);

                        if(pval>1 || pval<0) pval=1;

                        PvalLocus plocus;
                        plocus.raw_pval = pval;
                        plocus.coverage_factor=cover1;
                        plocus.meth_factor=countC1;
                        plocus.coverage_rest=cover2;
                        plocus.meth_rest=countC2;
                        plocus.chrom=chrom;
                        plocus.context=context_str[o1->context[j]];
                        plocus.position=o1->start[j];
                        plocus.sign=strand_str[o1->strand[j]];
                        if(ifp1->hdr->version & BM_ID)
                            plocus.name=o1->entryid[j];
                        else
                            plocus.name="";
                        plocus.pos = chrom_offset + 1 + plocus.position;
                        //meth_diff = fabs(meth_diff);
                        if(meth_diff>=methdiff && pval<=Pcutoff){
                    //#pragma omp critical
                            pvals.push_back(plocus);
                        }
                    //#pragma omp critical
                        if(meth_diff>=methdiff)
                        NcountC++;
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
}

int bm_overlap_mul_calp(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand, int Nsample1, Regression &full_regression, Regression &null_regression, vector<PvalLocus> &pvals, unsigned long &NcountC, \
    float Pcutoff, float methdiff, unsigned long chrom_offset, char *filtercontext){

//    bmOverlappingIntervals_t *o1;
    //fprintf(stderr, "process region %s %d %d\n", chrom, start, end);
    int slen = 1, i =0, j = 0;
    int* countM = (int *)malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0
    
    //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
    if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
    int total = 0, loci = 0;
    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) {
            fprintf(stderr, "Received an error somewhere!\n");
            return -1;
        };
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
                if(strcmp(filtercontext, "C")!=0){
                    if(ifp1s[i]->hdr->version & BM_CONTEXT){
                        if(strcmp(filtercontext, context_str[o1->context[j]])!=0){
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

    char **printmr = (char **)malloc(sizeof(char*)*(end-start+1));
    for(i=0;i<end-start+1;i++){
        printmr[i] = (char *)malloc(200+sizeifp*20);
    }
    char *tempchar = (char *)malloc(20);
    double pval = 0;

    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) {
            fprintf(stderr, "Received an error somewhere!\n");
            return -1;
        };
        if(o1->l){
            for(j=0; j<o1->l; j++) {
                loci = o1->start[j]-start;
                if(countM[loci] == sizeifp){
                    if(i == 0){ // file idx, file 0
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
                        sprintf(tempchar,"\t%u", o1->coverage[j]);
                        strcat(printmr[loci], tempchar);
                    }else{
                        sprintf(tempchar,"\t%f", o1->value[j]);
                        strcat(printmr[loci], tempchar);
                    }
                    
                    if(i==sizeifp-1){ // last i
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

    full_regression.props.chrom.clear();
    full_regression.props.chrom = chrom;
    string chromtmp = "";
    unsigned int countC1=0, cover1=0, countC2=0, cover2=0;
    int nsamples=0;
    for(i=0;i<end-start+1;i++){
        if(countM[i] == sizeifp){ // print merge results
            if(Nsample1 == 0) printf("%s\n",printmr[i]);
            else {
                //chr1    1958197 CG  -   6   8   6   8   0   5   0   5
                istringstream row_stream(printmr[i]);
                full_regression.props.position = 0;
                full_regression.props.strand.clear();
                full_regression.props.context.clear();
                full_regression.props.meth.clear();
                full_regression.props.total.clear();
                full_regression.props.name.clear();

                row_stream >> chromtmp;
                row_stream >> full_regression.props.position;
                row_stream >> full_regression.props.context;
                row_stream >> full_regression.props.strand;

                size_t total_count, meth_count;
                countC1=0, cover1=0, countC2=0, cover2=0;
                nsamples=0;
                while (row_stream >> meth_count >> total_count ) {
                    full_regression.props.total.push_back(total_count);
                    full_regression.props.meth.push_back(meth_count);
                    nsamples++;
                    if(nsamples <= Nsample1){
                        countC1 += meth_count;
                        cover1 += total_count;
                    }else{
                        countC2 += meth_count;
                        cover2 += total_count;
                    }
                }
                //printf("%d %d %d %d %d\n", i, meth_count, total_count, meth_count, total_count);
                if(geneid) { row_stream >> full_regression.props.name;}

                //cal p-value
                double meth_diff=fabs(double(countC1)/cover1 - double(countC2)/cover2);
                if(meth_diff>=methdiff) {
                    fit(full_regression);
                    null_regression.props = full_regression.props;
                    fit(null_regression);
                    pval = loglikratio_test(null_regression.max_loglik, full_regression.max_loglik);
                }else {
                    pval = 1;
                    NcountC++;
                    continue;
                }

                // If error occured in the fitting algorithm (i.e. p-val is nan or -nan).
                pval= ( (pval != pval) ? -1 : pval);
                //fprintf(stderr, "%s %d %s %f\n", chrom, full_regression.props.position, full_regression.props.context.c_str(), pval);
                if(meth_diff>=methdiff && pval<=Pcutoff){
                    PvalLocus plocus;
                    plocus.raw_pval = pval;
                    plocus.coverage_factor=cover1;
                    plocus.meth_factor=countC1;
                    plocus.coverage_rest=cover2;
                    plocus.meth_rest=countC2;
                    plocus.chrom=chrom;
                    plocus.context=full_regression.props.context;
                    plocus.position=full_regression.props.position;
                    plocus.sign=full_regression.props.strand;
                    if(geneid)
                        plocus.name=full_regression.props.name;
                    else
                        plocus.name="";
                    plocus.pos = chrom_offset + 1 + plocus.position;
                    pvals.push_back(plocus);
                }
                NcountC++;
            }
        }
    }

    //free(chrom);
    //bmDestroyOverlappingIntervals(o1);
    for(i=0;i<end-start+1;i++){
        free(printmr[i]);
    }
    free(tempchar);
    free(printmr);
    free(countM);
    return 0;
    
}

int bm_overlap_mul(binaMethFile_t **ifp1s, int sizeifp, char *chrom, int start, int end, uint8_t strand, int Nsample1){
    //fprintf(stderr, "process region %s %d %d\n", chrom, start, end);
    int slen = 1, i =0, j = 0;
    int* countM = (int *)malloc(sizeof(int)*(end-start+1));
    memset(countM, 0, sizeof(int)*(end-start+1)); // init 0

    //sscanf((const char *)regions[i], "%s:%d-%d", chrom, &start, &end);
    if(DEBUG>1) fprintf(stderr, "slen %d %d chrom %s %d %d %d", slen, i, chrom, start, end, slen);
    int total = 0, loci = 0;
    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) {
            fprintf(stderr, "Received an error somewhere!\n");
            return -1;
        };
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

    char **printmr = (char **)malloc(sizeof(char*)*(end-start+1));
    for(i=0;i<end-start+1;i++){
        printmr[i] = (char *)malloc(200+sizeifp*20);
    }
    char *tempchar = (char *)malloc(20);

    for(i=0;i<sizeifp;i++){
        bmOverlappingIntervals_t *o1;
        o1 = bmGetOverlappingIntervals(ifp1s[i], chrom, start, end+1);
        if(!o1) {
            fprintf(stderr, "Received an error somewhere!\n");
            return -1;
        };
        if(o1->l){
            for(j=0; j<o1->l; j++) {
                loci = o1->start[j]-start;
                if(countM[loci] == sizeifp){
                    if(i == 0){ // file idx
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
                        sprintf(tempchar,"\t%u", o1->coverage[j]);
                        strcat(printmr[loci], tempchar);
                    }else{
                        sprintf(tempchar,"\t%f", o1->value[j]);
                        strcat(printmr[loci], tempchar);
                    }

                    if(i==sizeifp-1){ // last i
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
        if(countM[i] == sizeifp){ // print merge results
            if(Nsample1 == 0) printf("%s\n",printmr[i]);
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


static inline double log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
};

// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
          gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
          (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
};


static double fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
};

static bool lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool lt_locus_raw_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.raw_pval < r2.raw_pval;
}

static bool ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.pos < r2.pos;
}
////BH
void adjust(vector<PvalLocus> &loci, unsigned long NcountC) {

      std::sort(loci.begin(), loci.end(), lt_locus_raw_pval);

      for (size_t ind = 0; ind < loci.size(); ++ind) 
      {
        const double current_score = loci[ind].raw_pval;
        //Assign a new one.
        //const double adjust_pval 
        loci[ind].adjust_pval = loci.size()*current_score/(ind + 1); 
        //loci[ind].adjust_pval = NcountC*current_score/(ind + 1 + NcountC - loci.size()); 
        //loci[ind].adjust_pval = current_score*NcountC/(ind + 1);
        //= adjust_pval;
      }

      for (vector<PvalLocus>::reverse_iterator it = loci.rbegin() + 1; it != loci.rend(); ++it) 
      {
        const PvalLocus &prev_locus = *(it - 1);
        PvalLocus &cur_locus = *(it);//cout << cur_locus.adjust_pval <<endl;
        cur_locus.adjust_pval = std::min(prev_locus.adjust_pval, cur_locus.adjust_pval);
      }

      for (vector<PvalLocus>::iterator it = loci.begin(); it != loci.end(); ++it) {
        PvalLocus &cur_locus = *(it); 
        if (cur_locus.adjust_pval > 1.0)
          cur_locus.adjust_pval = 1.0;
      }
	
      // Restore original order
      std::sort(loci.begin(), loci.end(), ls_locus_position); 
}

void Print_dm_result(const vector<PvalLocus> &pval_loci, ostream &output_encoding, double fdr_cutoff, double Pcutoff, double methdiff, string dmr_outfile, 
    bool singleAuto, int mindmc, int mindmcdis, int maxdmrlen, bool geneid) {
  string record, chrom, context, sign;string name;
  int position, coverage_factor, meth_factor, coverage_rest, meth_rest; //size_t
  double pval;double adjust_pvalue;
  ofstream OutFiledmr;
  if(singleAuto) OutFiledmr.open(dmr_outfile.c_str());

  vector<PvalLocus>::const_iterator cur_locus_iter = pval_loci.begin();
  int Ndmc = 0; int prevdmc = 0; int dmrtotalLen = 0; string prevchrom = "null";
  int dmrC1 = 0; int dmrCover1 = 0; int dmrC2 = 0; int dmrCover2 = 0;
  int dmrstart = 0; int dmrend = 0;
  int hyperNdmc=0,hypoNdmc=0;
  for(unsigned i=0;i<pval_loci.size();) 
  {
    pval= cur_locus_iter->raw_pval;
    adjust_pvalue=cur_locus_iter->adjust_pval; 
    coverage_factor=cur_locus_iter->coverage_factor;
    meth_factor=cur_locus_iter->meth_factor;
    coverage_rest=cur_locus_iter->coverage_rest;
    meth_rest=cur_locus_iter->meth_rest;
    position = cur_locus_iter->position;
    chrom = cur_locus_iter->chrom;
    sign = cur_locus_iter-> sign ; context = cur_locus_iter->context ; //std::string state="Faild";
    std::string name=cur_locus_iter->name.c_str();
    double meth_diff=fabs(double(meth_factor)/coverage_factor - double(meth_rest)/coverage_rest);
    double signed_methdiff=double(meth_factor)/coverage_factor - double(meth_rest)/coverage_rest;
    //if(!(adjust_pvalue < cutoff && (pval < Pcutoff || (pval < Pcutoff+0.05 && coverage_factor+coverage_rest<=50) ) && meth_diff >= methdiff ) ) {
    //fprintf(stderr, "111\n");
    if(!(adjust_pvalue <= fdr_cutoff && pval <= Pcutoff && meth_diff >= methdiff ) ) {
        cur_locus_iter++;
        i++;
        continue;
    };
    //fprintf(stderr, "222\n");
    
    output_encoding << chrom << "\t" << position << "\t" << sign << "\t"
                    << context << "\t" << pval << "\t";
    if (0 <= pval && pval <= 1) {
      output_encoding << adjust_pvalue << "\t";
      cur_locus_iter++;i++;
    } else {
      output_encoding << pval << "\t";
      cur_locus_iter++;i++;
    }
    
    //if(geneid) name=cur_locus_iter->name.c_str();
    if(geneid) output_encoding << meth_factor  << "\t" << coverage_factor << "\t"
                    << meth_rest  << "\t" << coverage_rest << "\t" << signed_methdiff << "\t" <<  name.c_str() << "\n" ;
    else {
        output_encoding << meth_factor  << "\t" << coverage_factor << "\t"<< meth_rest  << "\t" << coverage_rest << "\t" << signed_methdiff
    			<< "\n" ;
        if(singleAuto){
            if(prevchrom=="null") prevchrom = chrom;
//fprintf(stderr, "\nHHHH %d %d %d %d %d %d %d", Ndmc, mindmc, dmrstart, dmrend, prevdmc, position, position - prevdmc );
            if(dmrstart == 0){ 
                dmrstart = position;
                prevdmc = position;
                Ndmc=1 ;
                if(signed_methdiff > 0) {
                     hyperNdmc=1;
                     hypoNdmc=0;
                }
                else{
                     hypoNdmc=1;
                     hyperNdmc=0;
                }
            }
            else if(prevchrom != chrom) {
//fprintf(stderr, "\nzzzTTT %d %d %d %s %s", dmrtotalLen, dmrtotalLen + position - prevdmc, maxdmrlen, prevchrom.c_str(), chrom.c_str());
                if(Ndmc >= mindmc){
                    OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                    << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                    << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                }
                Ndmc = 0; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                dmrC1 = 0; dmrCover1 = 0; dmrC2 = 0; dmrCover2 = 0;
                dmrstart = position;
                prevchrom = chrom;
                hyperNdmc = 0; hypoNdmc = 0;
            }
            else if(position - prevdmc <= mindmcdis) {
        //fprintf(stderr, "\nTTT %d %d %d", dmrtotalLen, dmrtotalLen + position - prevdmc, maxdmrlen);
               if(dmrtotalLen + position - prevdmc <= maxdmrlen){ //未超出最大长度
                   Ndmc ++ ;
                   if(signed_methdiff > 0) hyperNdmc++;
                   else hypoNdmc++;
                   dmrend = position;
                   dmrtotalLen = dmrtotalLen + position - prevdmc;
                   prevdmc = position;
                   dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                   dmrC2 += meth_rest; dmrCover2 += coverage_rest;
               }else{ //超出了最大长度
                   if(Ndmc >= mindmc){
                       OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                        << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                        << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                   }
                   Ndmc = 1; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                   if(signed_methdiff > 0) {
                        hyperNdmc=1;
                        hypoNdmc=0;
                   }
                   else{
                        hypoNdmc=1;
                        hyperNdmc=0;
                   }
                  dmrC1 = 0; dmrCover1 = 0; dmrC2 = 0; dmrCover2 = 0;
                   dmrstart = position;
                   dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                   dmrC2 += meth_rest; dmrCover2 += coverage_rest;
               }
            }else{ // jian ge tai da, should be here
                if(Ndmc >= mindmc){
                     OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                     << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                     << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                }
                Ndmc = 1; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                if(signed_methdiff > 0) {
                     hyperNdmc=1;
                     hypoNdmc=0;
                }
                else{
                     hypoNdmc=1;
                     hyperNdmc=0;
                }
                dmrC1 = 0; dmrCover1 = 0; dmrC2 = 0; dmrCover2 = 0;
                dmrstart = position;
                dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                dmrC2 += meth_rest; dmrCover2 += coverage_rest;
            }
        }// end single auto
    } //end print
  } //end for
//printf("\n%d %d\n", Ndmc, mindmc);
  if(Ndmc >= mindmc){
     OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
     << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
     << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
     dmrC1 = 0; dmrCover1 = 0; dmrC2 = 0; dmrCover2 = 0;
  }
  return;
};

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
double loglikratio_test(double null_loglik, double full_loglik) {

  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2*(null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full model.
  // Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double pval = 1.0 - chisq_p;

  return pval;
}
