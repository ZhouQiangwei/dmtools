#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include "limits.h"
#include <map>
#include <algorithm>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include "binaMeth.h"
#include <stdlib.h>

#include "htslib/htslib/sam.h"
#include "htslib/htslib/hts.h"
#include <zlib.h>
#include "htslib/htslib/faidx.h"

#define MAX_HITS_ALLOWED 1500
#define CHROMSIZE 100
#define BATBUF 50000
#define MAXTAG 500

#define QLIMIT_FLOAT 30.0f
//#define max(a,b) ( ((a)>(b)) ? (a):(b) )
//#define min(a,b) ( ((a)>(b)) ? (b):(a) )

const float QLIMIT=QLIMIT_FLOAT;

struct Mismatch
{
	char Base;
	int Pos;
};
struct Gene_Hash
{
	char* Genome;
	int Index;
};
struct Methy_Hash
{
	int *plusMethylated,*plusUnMethylated;//,*plusCover
	int *plusG,*plusA;
	int *NegMethylated,*NegUnMethylated;//,*NegCover
	int *NegG,*NegA;
	//int *MethContext;
	int Index;
};
struct Offset_Record
{
	char Genome[1000];
	unsigned Offset;
} Temp_OR; 
typedef struct {
   Gene_Hash* Genome_List;
   Offset_Record* Genome_Offsets;
   char* Org_Genome;
   char* Marked_Genome;
   char* Marked_GenomeE;
   Methy_Hash Methy_List;
   FILE* OUTFILE;
   FILE* samINFILE;
   samFile* BamInFile;
   bam1_t *b;
   //bam_header_t *header;
   bam_hdr_t *header;
   int ThreadID;
   off64_t File_Size;
   char* processChr;
   char* INbamfilename;
} ARGS;
bool RELESEM = false;
bool printheader = true;
using namespace std;
bool Collision=false;
map <string,int> String_Hash;
float ENTROPY_CUTOFF=0;
///g++ ./src/split.cpp -o ./src/split -lpthread
//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
int fprintf_time(FILE *stream, const char *format, ...);
int strSearch(char *str1,char *str2);
int Get_ED(std::string & S);
char* replace(char*src, char*sub, char*dst);
void MinandSec(unsigned a[MAX_HITS_ALLOWED][10], int left, int right, int&min, int&second,int & loci);
//void MaxandSec(int a[], int left, int right, int & max, int & second);
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
void fetch(const char *str, char c1, char c2, char *buf);
string&   replace_all(string&   str,const   string&   old_value,const   string&   new_value);
std::string m_replace(std::string str,std::string pattern,std::string dstPattern,int count=-1);
char *fastStrcat(char *s, char *t);
float Pr(float Q);
void Build_Pow10();
void initMem(char* temp,int size,char* temp2);
void *Process_read(void *arg);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%  ++++++++]\n");//progress bar....
#define random(x) (rand()%x)
inline unsigned char Hash(char* S);
int Get_String_Length(FILE* INFILE);
string remove_soft_split(string & cigar,int & Read_Len,int & pos);

void ReverseC_Context(char* Dest,const char* seq,int & stringlength);
//---------------------------------------------------------------------------------------------------------------
unsigned Total_Reads_all;
unsigned long Total_Reads=0, Total_mapped = 0, forward_mapped = 0, reverse_mapped = 0;
//----mapping count
unsigned Tot_Unique_Org=0;//total unique hits obtained
unsigned ALL_MAP_Org=0;
unsigned Tot_Unique_Remdup=0;//total unique hits obtained after removing dups...
unsigned ALL_Map_Remdup=0;
float UPPER_MAX_MISMATCH=0.1;
bool REMOVE_DUP=false; //true; //true to removeDup, false will not remove PCR-dup
unsigned Mismatch_Qual[255][255][255]; //[readLength][255][255]
int QualCut=20;
const int POWLIMIT=300;
float POW10[POWLIMIT];
int QUALITYCONVERSIONFACTOR=33;
//----------
bool Methratio=false;
//-------meth count
unsigned long non_met_CG=0;
unsigned long met_CG=0;
unsigned long non_met_CHG=0;
unsigned long met_CHG=0;
unsigned long non_met_CHH=0;
unsigned long met_CHH=0;
bool bamformat=false;

int Sam=1;//1 true 0 false
unsigned Number_of_Tags = 0;
//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
char Char2Comp[255];
const int INITIAL_PROGRESS_READS =10000;//progress bar initial count..
float calc_Entropy(string readString, int L);
void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L);
void print_meth_tofile(int genome_id, ARGS* args);
bool SamSeqBeforeBS=false;
int RegionBins=1000;
long longestChr = 0;
string processingchr = "NULL";
string newchr = "NULL";
FILE* METHOUTFILE;
//FILE* BINsOUTFILE;
//FILE* REGION_OUT_CG;
//FILE* REGION_OUT_CHG;
//FILE* REGION_OUT_CHH;
int printtxt=0;
//FILE* LOC_OUT_CG;
//FILE* LOC_OUT_CHG;
//FILE* LOC_OUT_CHH;
int Mcoverage=4;
int maxcoverage=500;
int binspan=50000;
int nCs=1;
//C coverage
unsigned long totalC = 0;
unsigned long totalG = 0;
unsigned long CoverC[16]={0};
//
unsigned long M=0,Mh=0,H_AllC=0,hU=0,U=0;
unsigned long M_CG=0,Mh_CG=0,H_CG=0,hU_CG=0,U_CG=0;
unsigned long mCGdensity[100]={0},mCHGdensity[100]={0},mCHHdensity[100]={0};
unsigned long plus_mCGcount=0,plus_mCHGcount=0,plus_mCHHcount=0;
unsigned long plusCGcount=0,plusCHGcount=0,plusCHHcount=0;
unsigned long Neg_mCGcount=0,Neg_mCHGcount=0,Neg_mCHHcount=0;
unsigned long NegCGcount=0,NegCHGcount=0,NegCHHcount=0;
string Prefix="None";
#define MAX_LINE_PRINT 1000000
#define BUFSIZE 1000000
#define MAX_CHROM 1000000
binaMethFile_t *fp;
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	
	const char* Help_String="Command Format :  dmtools bam2dm [options] -g genome.fa -i/-b <SamfileSorted/BamfileSorted> -m <methratio dm outfile prefix>\n"
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
		"\t-nC		             >= <INT> nCs per region. default:1\n"
		"\t-r|--remove_dup       REMOVE_DUP, default:false\n"
        "\t--mrtxt               print prefix.methratio.txt file\n"
        "\t [DM format] paramaters\n"
        "\t-C                    print coverage\n"
        "\t-S                    print strand\n"
        "\t--Cx                  print context\n"
        "\t-E                    print end\n"
        "\t--Id                  print ID\n"
        "\t--zl                  The maximum number of zoom levels. [1-10], default: 2\n"
        "\t-i|--input            Sam format file, sorted by chrom.\n"
        "\t-h|--help";
/*
		"\t-as [0/1]             If print calculated alignment reads in sam/bam file. default:0\n"
		"\t-f|--sam              f for sam format outfile contain methState. default: sam format.\n"
		"\t--sam-seq-beforeBS    Converting BS read to the genome sequences.\n"
		"\t-h|--help";
*/

	Char2Comp['A']=Char2Comp['a']='T';
	Char2Comp['C']=Char2Comp['c']='G';
	Char2Comp['G']=Char2Comp['g']='C';
	Char2Comp['T']=Char2Comp['t']='A';
	Char2Comp['N']=Char2Comp['n']='N';
	int Genome_CountX=0;
	char Output_Name[100];
	strcpy(Output_Name, "None");
	//string Prefix="None";
	string methOutfileName;
	string Geno;
//	int Current_Option=0;
	int InFileStart=0,InFileEnd=0;
	char *InFile;
	string binsOutfileName="";
	
//	int Par_Count=0;
	std::string CMD;
	int AlignSum = 0;
	uint32_t write_type = 0x8000; //0xf1ff
	int chrlenf = -1;
	string chrsizefile;
	int zoomlevel = 2;
    char processChr[1000];
    strcpy(processChr, "NAN-mm");

	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-f") ||!strcmp(argv[i], "--sam")  )
		{
			strcpy(Output_Name, argv[++i]);
			Sam=1;
		}
		else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "--chrom"))
        {
            strcpy(processChr, argv[++i]);
        }else if(!strcmp(argv[i], "-Q"))
		{
			QualCut=atoi(argv[++i]);
        }
		else if(!strcmp(argv[i], "--sam-seq-beforeBS"))
		{
			SamSeqBeforeBS=true;
		}
		else if(!strcmp(argv[i], "-r") || !strcmp(argv[i], "--remove_dup"))
		{
			REMOVE_DUP=true;
		}else if(!strcmp(argv[i], "-n") || !strcmp(argv[i], "--Nmismatch"))
		{    
			UPPER_MAX_MISMATCH=atof(argv[++i]);
			if(UPPER_MAX_MISMATCH>1){
				fprintf(stderr, "\nError defined mismatch paramater, should be 0-1.\n");
			}
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
        }
		else if(!strcmp(argv[i], "-R") || !strcmp(argv[i], "--Regions"))
		{
			RegionBins=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "--binsfile"))
		{
			binsOutfileName=argv[++i];
		}
		else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--step"))
		{
			binspan=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i],"-m") || !strcmp(argv[i],"--methratio"))
		{
			Methratio=true;
			Prefix=argv[++i];
		}
		else if(!strcmp(argv[i],"-c") || !strcmp(argv[i],"--coverage"))
		{
			Mcoverage=atoi(argv[++i]);
        }else if(!strcmp(argv[i],"--mrtxt"))
		{
			fprintf(stderr, "Print methratio.txt\n");
			printtxt=1;
        }else if(!strcmp(argv[i],"-as"))
		{
			AlignSum=atoi(argv[++i]);
        }else if(!strcmp(argv[i],"--maxcoverage"))
		{
			maxcoverage=atoi(argv[++i]);
        }else if(!strcmp(argv[i],"-nC"))
		{
			nCs=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-i") || !strcmp(argv[i], "--input"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-')
			{
				i++;
				continue;
			}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
		}else if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--binput"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-')
			{
				i++;
				continue;
			}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
			bamformat=true;
		}
		else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
			printf("\n%s\n",Help_String);
                        exit(0);
		}
		else
		{
			printf("\nError: %s\n\n%s \n",argv[i],Help_String);
			exit(0);
		}
		
	}
	for(int i = 0; i < argc; i++) {CMD.append(argv[i]); CMD.append(" ");}
	
	if(argc<=3) printf("%s \n",Help_String);
	if (argc >3  && InFileStart) 
	{
		printf("\n[DMtools::calmeth] v1.0\n");
		string log;
                log=Prefix;
                log+=".methlog.txt";
		FILE* OUTLOG=File_Open(log.c_str(),"w");
		string mCdensity;
                mCdensity=Prefix;  
                mCdensity+=".mCdensity.txt";

                string NCcoverage;
                NCcoverage=Prefix;
                NCcoverage+=".NCcoverage.txt";

		string mCcatero;
                mCcatero=Prefix;
                mCcatero+=".mCcatero.txt";

		printf("[DM::calmeth] Coverage and validC: %d %d, %d\n", Mcoverage, maxcoverage, nCs);	
		try
		{
			time(&Start_Time);
			Build_Pow10();

			
			char **chroms = (char **)malloc(sizeof(char*)*MAX_CHROM);
			//if(!chroms) goto error;
			uint32_t *chrLens = (uint32_t *)malloc(sizeof(uint32_t) * MAX_CHROM);
			
			//string G=Geno;G+=".bin";
			//string chrLenf=Geno;chrLenf+=".len"; //ann.location

			FILE* GenomeFILE=File_Open(Geno.c_str(),"r");
			printf("[DM::calmeth] Loading genome sequence : %s\n", Geno.c_str());
			ARGS args;
			fseek(GenomeFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(GenomeFILE);rewind(GenomeFILE);//load original genome..
			args.Org_Genome=new char[Genome_Size];if(!args.Org_Genome) throw("Insufficient memory to load genome..\n"); 
			char readBuffer[BUFSIZE];
			int Genome_Count=0, lines = 0, chrlen = 0;
			char* token;
			char seps[] = " \t\n\r";
			char chrName[1000]; int perlen = 0;
			Genome_Size = 0;
			args.Genome_Offsets = new Offset_Record[MAX_CHROM];
			args.Org_Genome[0] = '\0';
			char *genometmp = args.Org_Genome;int kk = 0;
			
			while (fgets(readBuffer,BUFSIZE,GenomeFILE)!=0)//count genomes..
			{
				if(strlen(readBuffer) >= BUFSIZE - 1) {
					fprintf(stderr, "Too many characters in one row! Try to split the long row into several short rows (fewer than %d characters per row).\n", BUFSIZE);
					exit(1);
				}
				if(readBuffer[0] == '>') {
					if(lines > 0) {
						fprintf(stderr, "[DM::calmeth] loaded and store %s %d %d\n", chrName, chrlen, Genome_Count);
						lines = 0;
                        if(strcmp(processChr, "NAN-mm") == 0) {
						    if(longestChr < chrlen) longestChr = chrlen;
                        }else if(strcmp(chrName, processChr) == 0) {
                            longestChr = chrlen;
                        }

						strcpy(args.Genome_Offsets[Genome_Count].Genome, chrName);
						args.Genome_Offsets[Genome_Count].Offset = chrlen;
						
						chrLens[Genome_Count] = chrlen;
						chroms[Genome_Count] = strdup(chrName);
					}

					// Save name
					token = strtok(readBuffer + 1, seps);	
					strcpy(chrName, token);

					chrlen = 0;
					Genome_Count++;
					if(Genome_Count>MAX_CHROM){
						fprintf(stderr, "\nchrom count bigger than 1000000, please -x-\n");
						exit(0);
					}
				}
				else if(strcmp(processChr, "NAN-mm") == 0 || strcmp(chrName, processChr) == 0) {
					// Substract \n
					perlen = strlen(readBuffer) - 1;
					chrlen += perlen;
					//strcat(args.Org_Genome, readBuffer);
					genometmp = fastStrcat(genometmp, readBuffer);
					Genome_Size+=perlen;
					//fprintf(stderr, "\n-- %d\n", Genome_Size);
					for(kk=0;kk<perlen;kk++){
						if(toupper(readBuffer[kk]) == 'G'){totalG++;}
						if(toupper(readBuffer[kk]) == 'C'){totalC++;}
					}
				}
				lines++;
			}
			if(lines > 0) {
                if(strcmp(processChr, "NAN-mm") == 0) {
                    if(longestChr < chrlen) longestChr = chrlen;
                }else if(strcmp(chrName, processChr) == 0) {
                    longestChr = chrlen;
                }

				strcpy(args.Genome_Offsets[Genome_Count].Genome, chrName);
				args.Genome_Offsets[Genome_Count].Offset = chrlen;
				
				chrLens[Genome_Count] = chrlen;
				chroms[Genome_Count] = strdup(chrName);
				//fprintf(stderr, "%s\t%d\n", chrName, chrlen);
			}

            if(longestChr == 0){
                fprintf(stderr, "Warining: undetected %s in genome file", processChr);
                exit(0);
            }

			fprintf(stderr, "[DM::calmeth] len count GCcount %lld %d %ld\n", Genome_Size, Genome_Count, totalC+totalG);
			//exit(0);
			//if(!fread(args.Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");
			if(REMOVE_DUP){
				args.Marked_Genome=new char[Genome_Size+1];if(!args.Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 
				args.Marked_GenomeE=new char[Genome_Size+1];if(!args.Marked_GenomeE) throw("Insufficient memory to Mark genome..\n"); 
			}
			
			args.OUTFILE = NULL;
			assert(longestChr>0);
			printf("[DM::calmeth] Longest chr: %d\n",longestChr);
			if(Sam && strcmp(Output_Name,"None") ) args.OUTFILE=File_Open(Output_Name,"w");

			char* Split_Point=args.Org_Genome;//split and write...
			args.Genome_List = new Gene_Hash[Genome_Count];
			//args.Methy_List = new Methy_Hash;
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				String_Hash[args.Genome_Offsets[i].Genome]=i;
				//unsigned char H=Hash(Genome_Offsets[i].Genome);
				if(i==0) args.Genome_List[i].Genome = (Split_Point+=0);
				else args.Genome_List[i].Genome=(Split_Point+=args.Genome_Offsets[i-1].Offset);
				args.Genome_List[i].Index=i;
				//meth ini
			}

			if(Methratio)
			{
				try{
				args.Methy_List.plusG = new int[longestChr];
				args.Methy_List.plusA = new int[longestChr];
				args.Methy_List.NegG = new int[longestChr];
				args.Methy_List.NegA = new int[longestChr];
				//args.Methy_List[i].MethContext =new int[args.Genome_Offsets[i+1].Offset];
				args.Methy_List.plusMethylated = new int[longestChr];
				args.Methy_List.plusUnMethylated = new int[longestChr];
				//Methy_List[i].plusCover = new int[Genome_Offsets[i+1].Offset];
				args.Methy_List.NegMethylated = new int[longestChr];
				args.Methy_List.NegUnMethylated = new int[longestChr];
				//Methy_List[i].NegCover = new int[Genome_Offsets[i+1].Offset];
				//=========Genome_List[i].Genome;
				}catch(std::bad_alloc){
					fprintf(stderr, "\nbad alloc in main array!\n");
				}
			}
			
			fclose(GenomeFILE);
			//fclose(Location_File);

			////read file
			
			if(Methratio){
				if(printtxt == 1){
					methOutfileName=Prefix;
                    methOutfileName+=".methratio.txt";
					METHOUTFILE=File_Open(methOutfileName.c_str(),"w");
					fprintf(METHOUTFILE,"#chromsome\tloci\tstrand\tcontext\tC_count\tCT_count\tmethRatio\teff_CT_count\trev_G_count\trev_GA_count\tMethContext\t5context\n");
				}

				fp = NULL;
				methOutfileName=Prefix;
                methOutfileName+=".methratio.dm";

				if(bmInit(1<<17) != 0) {
        		    fprintf(stderr, "Received an error in dmInit\n");
		            return 1;
        		}
				fp = (binaMethFile_t*)bmOpen((char*)methOutfileName.data(), NULL, "w");
	        	fp->type = write_type;
				if(!fp) {
					fprintf(stderr, "An error occurred while opening example_output.dm for writingn\n");
					return 1;
				}

				/*
				//bin 
				binsOutfileName=Prefix;binsOutfileName+=".methBins.txt";
				BINsOUTFILE=File_Open(binsOutfileName.c_str(),"w");
				//---DMR
				string RegionOutFiles_CG=Prefix;RegionOutFiles_CG+="_Region.CG.txt";
				string RegionOutFiles_CHG=Prefix;RegionOutFiles_CHG+="_Region.CHG.txt";
				string RegionOutFiles_CHH=Prefix;RegionOutFiles_CHH+="_Region.CHH.txt";
				REGION_OUT_CG=File_Open(RegionOutFiles_CG.c_str(),"w");//RegionBins RegionOutFiles
				REGION_OUT_CHG=File_Open(RegionOutFiles_CHG.c_str(),"w");
				REGION_OUT_CHH=File_Open(RegionOutFiles_CHH.c_str(),"w");
				//--------DMR---------------//
				//----------DMC
				string LociOutFiles_CG=Prefix;LociOutFiles_CG+="_loci.CG.txt";
				string LociOutFiles_CHG=Prefix;LociOutFiles_CHG+="_loci.CHG.txt";
				string LociOutFiles_CHH=Prefix;LociOutFiles_CHH+="_loci.CHH.txt";
				LOC_OUT_CG=File_Open(LociOutFiles_CG.c_str(),"w");//RegionBins RegionOutFiles
				LOC_OUT_CHG=File_Open(LociOutFiles_CHG.c_str(),"w");
				LOC_OUT_CHH=File_Open(LociOutFiles_CHH.c_str(),"w");
				*/
			}
                        args.INbamfilename = new char[1000];
			for(int f=InFileStart;f<=InFileEnd;f++)
			{
				printf("[DM::calmeth] Processing %d out of %d. File: %s, %d\n\n", f-InFileStart+1,InFileEnd-InFileStart+1, argv[f], bamformat);
				//fseek(args.INFILE, 0L, SEEK_END);args.File_Size=ftello64(args.INFILE);rewind(args.INFILE);
				char s2t[BATBUF];
				//if(args.OUTFILE!=NULL)
 				{ // && printheader
					//samfile_t *bamin = 0;
					args.BamInFile = 0;
					args.header;
                                        if(bamformat)
                                        {
                                             strcpy(args.INbamfilename, argv[f]);
						if ((args.BamInFile = sam_open(argv[f], "rb")) == 0) {
								fprintf(stderr, "fail to open \"%s\" for reading.\n", argv[f]);
						}
						args.header = sam_hdr_read((samFile*)args.BamInFile);
						if(InFileStart && args.OUTFILE!=NULL) fprintf(args.OUTFILE,"%s",args.header->text);
						int chrLenfile = 0;
						if(chrLenfile == 0 && args.header){
							char *chrom = (char*)malloc(100);
							unsigned long chrlen = 0;
							//fprintf(stderr, "len %d", args.header->n_targets);
							Genome_Count = args.header->n_targets;
							for(int k=0;k<args.header->n_targets;k++){
								//fprintf(stderr, "%s %d\n", args.header->target_name[k], args.header->target_len[k]);
								chrLens[k] = args.header->target_len[k];
								chroms[k] = strdup(args.header->target_name[k]);
								continue;
							}
						}else{
							fprintf(stderr, "Can not find header in bam file, please provide chrom size file with --cs paramater\nNote. The sequence of chromosome names shall be consistent with that in BAM file.\n");
							exit(0);
						}
                    }
                    else if(f==InFileStart){
                    	args.samINFILE=File_Open(argv[f],"r");
						Genome_Count = 0;
						while (fgets(s2t,BATBUF,args.samINFILE)!=0 ){
							char *chrom = (char*)malloc(100);
							unsigned long chrlen = 0;
							if(s2t[0]=='@')
							{    
									//s2t[BATBUF]='\0';s2t[BATBUF-1]='\n'; //这个会报错没搞清
									if(Sam && args.OUTFILE!=NULL){
										fprintf(args.OUTFILE,"%s",s2t);
									}
									int chrLenfile = 0;
									if(chrLenfile == 0){
										if(DEBUG>1) fprintf(stderr, "len %s %d\n", s2t, Genome_Count);
										sscanf(s2t, "@SQ\t%*[^:]:%s\t%*[^:]:%ld", chrom, &chrlen);
										if(DEBUG>1) fprintf(stderr, "lenxxxx %d\n", chrlen);
										if(chrlen != 0){
											if(DEBUG>1) fprintf(stderr, "len2 %s %d\n", chrom, chrlen);
											chrLens[Genome_Count] = chrlen;
											chroms[Genome_Count] = chrom;
											Genome_Count++;
										}
									}else{
										fprintf(stderr, "Can not find header in sam file, please provide chrom size file with --cs paramater\nNote. The sequence of chromosome names shall be consistent with that in SAM file.\n");
										exit(0);
									}
									continue;
							}else
								break;
							rewind(args.samINFILE);
						}
					}
				}
				if(f==InFileStart){
					//Allow up to 10 zoom levels, though fewer will be used in practice
					if(bmCreateHdr(fp, zoomlevel)) exit(0);
					//Create the chromosome lists
					fp->cl = bmCreateChromList(chroms, chrLens, Genome_Count); //2
					if(!fp->cl) exit(0);
					//Write the header
					if(bmWriteHdr(fp)) exit(0);
					//Some example methlevel
					if(DEBUG>1) fprintf(stderr, "====HHH type %d\n", fp->type);
				}
				//nothreads
                args.processChr=new char[1000];
                strcpy(args.processChr, processChr);
				Process_read(&args);
				Done_Progress();
				if(!bamformat) fclose(args.samINFILE);
            	if(bamformat) 
            	{
                	bam_hdr_destroy(args.header);
                	sam_close(args.BamInFile);
            	}
			}

			map<string, int>::iterator iter;
			int H = -1;
			iter = String_Hash.find(newchr.c_str());
			if(iter != String_Hash.end()){
				fprintf(stderr, "Print output of %s\n", newchr.c_str());
				H = iter->second;
				print_meth_tofile(H, &args);
			}
			//
			 if(Methratio){
				if(printtxt == 1){
					fclose(METHOUTFILE);
				}
				/*
				fclose(BINsOUTFILE);
				fclose(REGION_OUT_CG);
				fclose(REGION_OUT_CHG);
				fclose(REGION_OUT_CHH);

				fclose(LOC_OUT_CG);
				fclose(LOC_OUT_CHG);
				fclose(LOC_OUT_CHH);
				*/
			 }
			//Print_Mismatch_Quality(OUTFILE_MM_QUALITY, Read_Len);
	//}
			printf("genome process done!\n");
			//if(args.OUTFILE!=NULL) fclose(args.OUTFILE);
			if(Sam && strcmp(Output_Name,"None") ) fclose(args.OUTFILE);
			
			printf("Raw count of Met_C in CG:\t%lu\n",met_CG);
			printf("Raw count of Non_Met_C in CG:\t%lu\n",non_met_CG);
			printf("Raw count of Met_C in CHG:\t%lu\n",met_CHG);
			printf("Raw count of Non_Met_C in CHG:\t%lu\n",non_met_CHG);
			printf("Raw count of Met_C in CHH:\t%lu\n",met_CHH);
			printf("Raw count of Non_Met_C in CHH:\t%lu\n",non_met_CHH);
			fprintf(OUTLOG,"Case\tValue\n");
			fprintf(OUTLOG,"Raw count of Met_C in CG:\t%lu\n",met_CG);
			fprintf(OUTLOG,"Raw count of Non_Met_C in CG:\t%lu\n",non_met_CG);
			fprintf(OUTLOG,"Raw count of Met_C in CHG:\t%lu\n",met_CHG);
			fprintf(OUTLOG,"Raw count of Non_Met_C in CHG:\t%lu\n",non_met_CHG);
			fprintf(OUTLOG,"Raw count of Met_C in CHH:\t%lu\n",met_CHH);
			fprintf(OUTLOG,"Raw count of Non_Met_C in CHH:\t%lu\n",non_met_CHH);
			printf("[CpG]\tM: %lu Mh: %lu H: %lu hU: %lu U: %lu\n",M_CG,Mh_CG,H_CG,hU_CG,U_CG);
            printf("\n[mC]\tM: %lu Mh: %lu H: %lu hU: %lu U: %lu\n",M,Mh,H_AllC,hU,U);
            fprintf(OUTLOG,"[CpG]\tM: %lu Mh: %lu H: %lu hU: %lu U: %lu\n",M_CG,Mh_CG,H_CG,hU_CG,U_CG);
            fprintf(OUTLOG,"[mC]\tM: %lu Mh: %lu H: %lu hU: %lu U: %lu\n",M,Mh,H_AllC,hU,U);
			if(Methratio)
			{
				FILE* fp_NCcoverage=File_Open(NCcoverage.c_str(),"w");
				fprintf(fp_NCcoverage,"NCcovergae");
				for(int i=0;i<15;i++)
				{
						fprintf(fp_NCcoverage,"\t%lu",CoverC[i]);
				}
				fprintf(fp_NCcoverage,"\nPercent_of_covergae");
				for(int i=0;i<15;i++)
				{
						fprintf(fp_NCcoverage,"\t%f",((double)CoverC[i])/(totalC+totalG));
				}
				fclose(fp_NCcoverage);

				FILE* mC_DENSITY=File_Open(mCdensity.c_str(),"w");
				fprintf(mC_DENSITY,"mCG");
				for(int i=0;i<100;i++)
				{
					fprintf(mC_DENSITY,"\t%lu",mCGdensity[i]);
				}
				fprintf(mC_DENSITY,"\nmCHG");
				for(int i=0;i<100;i++)
				{
					fprintf(mC_DENSITY,"\t%lu",mCHGdensity[i]);
				}
				fprintf(mC_DENSITY,"\nmCHH");
				for(int i=0;i<100;i++)
				{
					fprintf(mC_DENSITY,"\t%lu",mCHHdensity[i]);
				}
				fclose(mC_DENSITY);

				FILE* mC_catero=File_Open(mCcatero.c_str(),"w");
				fprintf(mC_catero,"M\t%lu\nMh\t%lu\nH\t%lu\nhU\t%lu\nU\t%lu",M,Mh,H_AllC,hU,U);
				fprintf(mC_catero,"\nCpG_M\t%lu\nCpG_Mh\t%lu\nCpG_H\t%lu\nCpG_hU\t%lu\nCpG_U\t%lu",M_CG,Mh_CG,H_CG,hU_CG,U_CG);
				fprintf(mC_catero,"\nmC\t%f\n", double (1*(plus_mCGcount+plus_mCHGcount+plus_mCHHcount+Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount))/double (plusCGcount+plusCHGcount+plusCHHcount+NegCGcount+NegCHGcount+NegCHHcount) );
				fprintf(mC_catero,"mCG\t%f\n", double (1*(plus_mCGcount+Neg_mCGcount))/double (plusCGcount+NegCGcount));
				fprintf(mC_catero,"mCHG\t%f\n", double (1*(plus_mCHGcount+Neg_mCHGcount))/double (plusCHGcount+NegCHGcount));
				fprintf(mC_catero,"mCHH\t%f\n", double (1*(plus_mCHHcount+Neg_mCHHcount))/double (plusCHHcount+NegCHHcount) );

				fclose(mC_catero);
				//+
				printf("\nStrand+ :\nmC/(C+T) {%ld / %ld} = %f% \n",(plus_mCGcount+plus_mCHGcount+plus_mCHHcount),(plusCGcount+plusCHGcount+plusCHHcount),double (100*(plus_mCGcount+plus_mCHGcount+plus_mCHHcount))/(plusCGcount+plusCHGcount+plusCHHcount) );
				printf("mCG/(CG+TG) {%ld / %ld} = %f% \n",plus_mCGcount,plusCGcount,double (100*(plus_mCGcount))/(plusCGcount) );
				printf("mCHG/(CHG+THG) {%ld / %ld} = %f% \n",plus_mCHGcount,plusCHGcount,double (100*(plus_mCHGcount))/(plusCHGcount) );
				printf("mCHH/(CHH+THH) {%ld / %ld} = %f% \n",plus_mCHHcount,plusCHHcount,double (100*(plus_mCHHcount))/(plusCHHcount) );
				//-
				printf("\nStrand- :\nmC/(C+T) {%ld / %ld} = %f% \n",(Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount),(NegCGcount+NegCHGcount+NegCHHcount),double (100*(Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount))/double (NegCGcount+NegCHGcount+NegCHHcount) );
				printf("mCG/(CG+TG) {%ld / %ld} = %f% \n",Neg_mCGcount,NegCGcount,double (100*(Neg_mCGcount))/double (NegCGcount) );
				printf("mCHG/(CHG+THG) {%ld / %ld} = %f% \n",Neg_mCHGcount,NegCHGcount,double (100*(Neg_mCHGcount))/double (NegCHGcount) );
				printf("mCHH/(CHH+THH) {%ld / %ld} = %f% \n",Neg_mCHHcount,NegCHHcount,double (100*(Neg_mCHHcount))/double (NegCHHcount) );
				//+ -
				printf("\nStrand+- :\nmC/(C+T) {%ld / %ld} = %f% \n",(plus_mCGcount+plus_mCHGcount+plus_mCHHcount+Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount),(plusCGcount+plusCHGcount+plusCHHcount+NegCGcount+NegCHGcount+NegCHHcount),double (100*(plus_mCGcount+plus_mCHGcount+plus_mCHHcount+Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount))/double (plusCGcount+plusCHGcount+plusCHHcount+NegCGcount+NegCHGcount+NegCHHcount) );
				printf("mCG/(CG+TG) {%ld / %ld} = %f% \n",plus_mCGcount+Neg_mCGcount,plusCGcount+NegCGcount,double (100*(plus_mCGcount+Neg_mCGcount))/double (plusCGcount+NegCGcount) );
				printf("mCHG/(CHG+THG) {%ld / %ld} = %f% \n",plus_mCHGcount+Neg_mCHGcount,plusCHGcount+NegCHGcount,double (100*(plus_mCHGcount+Neg_mCHGcount))/double (plusCHGcount+NegCHGcount) );
				printf("mCHH/(CHH+THH) {%ld / %ld} = %f% \n",plus_mCHHcount+Neg_mCHHcount,plusCHHcount+NegCHHcount,double (100*(plus_mCHHcount+Neg_mCHHcount))/double (plusCHHcount+NegCHHcount) );
				fprintf(OUTLOG,"##%ld\t%ld\t%ld\n",totalC, totalG, totalC+totalG);
				//+
                fprintf(OUTLOG,"Strand\t+\nmC/(C+T)\t{%ld / %ld} = %f% \n",(plus_mCGcount+plus_mCHGcount+plus_mCHHcount),(plusCGcount+plusCHGcount+plusCHHcount),double (100*(plus_mCGcount+plus_mCHGcount+plus_mCHHcount))/(plusCGcount+plusCHGcount+plusCHHcount) );
                fprintf(OUTLOG,"mCG/(CG+TG)\t{%ld / %ld} = %f% \n",plus_mCGcount,plusCGcount,double (100*(plus_mCGcount))/(plusCGcount) );
                fprintf(OUTLOG,"mCHG/(CHG+THG)\t{%ld / %ld} = %f% \n",plus_mCHGcount,plusCHGcount,double (100*(plus_mCHGcount))/(plusCHGcount) );
                fprintf(OUTLOG,"mCHH/(CHH+THH)\t{%ld / %ld} = %f% \n",plus_mCHHcount,plusCHHcount,double (100*(plus_mCHHcount))/(plusCHHcount) );
                //-
                fprintf(OUTLOG,"Strand\t-\nmC/(C+T)\t{%ld / %ld} = %f% \n",(Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount),(NegCGcount+NegCHGcount+NegCHHcount),double (100*(Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount))/double (NegCGcount+NegCHGcount+NegCHHcount) );
                fprintf(OUTLOG,"mCG/(CG+TG)\t{%ld / %ld} = %f% \n",Neg_mCGcount,NegCGcount,double (100*(Neg_mCGcount))/double (NegCGcount) );
                fprintf(OUTLOG,"mCHG/(CHG+THG)\t{%ld / %ld} = %f% \n",Neg_mCHGcount,NegCHGcount,double (100*(Neg_mCHGcount))/double (NegCHGcount) );
                fprintf(OUTLOG,"mCHH/(CHH+THH)\t{%ld / %ld} = %f% \n",Neg_mCHHcount,NegCHHcount,double (100*(Neg_mCHHcount))/double (NegCHHcount) );
				//+ -
				fprintf(OUTLOG,"Strand\t+-\nmC/(C+T)\t{%ld / %ld} = %f% \n",(plus_mCGcount+plus_mCHGcount+plus_mCHHcount+Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount),(plusCGcount+plusCHGcount+plusCHHcount+NegCGcount+NegCHGcount+NegCHHcount),double (100*(plus_mCGcount+plus_mCHGcount+plus_mCHHcount+Neg_mCGcount+Neg_mCHGcount+Neg_mCHHcount))/double (plusCGcount+plusCHGcount+plusCHHcount+NegCGcount+NegCHGcount+NegCHHcount) );
				fprintf(OUTLOG,"mCG/(CG+TG)\t{%ld / %ld} = %f% \n",plus_mCGcount+Neg_mCGcount,plusCGcount+NegCGcount,double (100*(plus_mCGcount+Neg_mCGcount))/double (plusCGcount+NegCGcount));
				fprintf(OUTLOG,"mCHG/(CHG+THG)\t{%ld / %ld} = %f% \n",plus_mCHGcount+Neg_mCHGcount,plusCHGcount+NegCHGcount,double (100*(plus_mCHGcount+Neg_mCHGcount))/double (plusCHGcount+NegCHGcount));
				fprintf(OUTLOG,"mCHH/(CHH+THH)\t{%ld / %ld} = %f% \n",plus_mCHHcount+Neg_mCHHcount,plusCHHcount+NegCHHcount,double (100*(plus_mCHHcount+Neg_mCHHcount))/double (plusCHHcount+NegCHHcount) );
			}

            // print align summary
			if(AlignSum){
				log=Prefix + ".align.log.txt";
				FILE* ALIGNLOG=File_Open(log.c_str(),"w");
				fprintf(ALIGNLOG, "Total_reads\t%u\n", Total_Reads);
				fprintf(ALIGNLOG, "Total_Mapped_Q%d\t%u\n", QualCut, Total_mapped);
				fprintf(ALIGNLOG, "Mapped_forword\t%u\n", forward_mapped);
				fprintf(ALIGNLOG, "Mapped_reverse\t%u\n", reverse_mapped);
				fclose(ALIGNLOG);
			}
			fprintf(stderr, "[DM::calmeth] dm closing\n");
        	bmClose(fp);
    	    fprintf(stderr, "[DM::calmeth] dm closed\n");
	        bmCleanup();
//delete
			fprintf(stderr, "[DM::calmeth] Done and release memory!\n");
			for(int i =0; i < MAX_CHROM; i++){
        	    free(chroms[i]);
    	    }
			free(chroms);free(chrLens);
			if(RELESEM){
				for ( int i=0;i<1;i++)
				{
					if(Methratio)
					{
						delete[] args.Methy_List.plusG;
						delete[] args.Methy_List.plusA;
						delete[] args.Methy_List.NegG;
						delete[] args.Methy_List.NegA;
						delete[] args.Methy_List.plusMethylated;
						delete[] args.Methy_List.plusUnMethylated;
						delete[] args.Methy_List.NegMethylated;
						delete[] args.Methy_List.NegUnMethylated;
					}
				}
				delete [] args.Genome_List;
				delete [] args.Genome_Offsets;
				delete [] args.Org_Genome;
			}
		}
		catch(char* Err)
		{
			printf(Err);
			fprintf(stderr, "\nError cigar\n");
			exit(-1);
		}

		time(&End_Time);printf("[DM::calmeth] Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	
	fclose(OUTLOG);
	}
	
}

int Read_Tag(FILE *INFILE,char s2t[],string hits[],int& cntHit)
{
	flockfile(INFILE);
	char Buf[BATBUF];
	if (!feof(INFILE))
	{
		Total_Reads_all++;
		fgets(s2t,BATBUF,INFILE);//read description..
		cntHit=0;
		fgets(Buf,BATBUF,INFILE);
		while(!feof(INFILE) && Buf[0]!='@')//while in the record.. 
		{
			hits[cntHit++]=Buf; //read hit info..
			fgets(Buf,BATBUF,INFILE);
			if(cntHit>=MAX_HITS_ALLOWED) 
			{
				cntHit=0;
		            while(!feof(INFILE) && Buf[0]!='@')
		            {
		                    fgets(Buf,BATBUF,INFILE);
		            }
				break;
			}
			assert(MAX_HITS_ALLOWED > cntHit);
		}
		funlockfile(INFILE);
		return 1;
	}
	else 
	{
		funlockfile(INFILE);
		return 0;
	}
}

char* process_cigar(const char* cig,int Read_Len)
{
	char temp[8];
	char* cigar_rm = new char[1000]();*cigar_rm=0;unsigned n=0;
	char* buffer_cigar=new char[100]();*buffer_cigar=0;
	while(*cig!='\0')
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				sprintf(buffer_cigar,"%dS",length);
				strcat(cigar_rm, buffer_cigar);
			}
            cig++;n=0;
		}else if(*cig=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				sprintf(buffer_cigar,"%dM",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
			char buf[1024]="\0";
			sprintf( buf , "%dM",Read_Len);
			char buf_tmp[1024]="\0";
			sprintf( buf_tmp, "%dM",length);
			if( !strcmp(buf,  buf_tmp ) )
			{
				if(!strcmp(buf, cig))
					break;
				else
				{
					strcpy(cigar_rm, buffer_cigar);
				}
				
			}
		}else if(*cig=='I')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1I");
				else if(length==2) sprintf(buffer_cigar,"2I");
				else sprintf(buffer_cigar,"%dI",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
		}else if(*cig=='D')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1D");
				else if(length==2) sprintf(buffer_cigar,"2D");
				else sprintf(buffer_cigar,"%dD",length);
				strcat(cigar_rm, buffer_cigar);
			}
			cig++;n=0;
		}else
		{
			printf(" --%d%c-- ",atoi(temp),*cig);
    		continue;//break;
		}
	}
	delete []cigar_rm;
	delete []buffer_cigar;
	return cigar_rm;
}
string getstring(char* seq, int l, int len){
	char tmp[10];
	for(int i=0; i<len; i++){
//printf("\ns %d %d %c\n", l, i, seq[l+i]);
		tmp[i] = seq[l+i];
	}
	if(len>0) tmp[len]='\0';
	return tmp;
}
//g++ bmtools/libBigWig.a bmtools/libBigWig.so calmeth.cpp -o calmeth -m64 -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz
char **chromsUse;
char **entryid;
uint32_t *starts;
uint32_t *pends;
float *values;
uint16_t *coverages;
uint8_t *strands;
uint8_t *contexts;
void print_meth_tofile(int genome_id, ARGS* args){
	if(Methratio)
	{
		fprintf(stderr, "[DM::calmeth] Start process chrom %d\n", genome_id);
		chromsUse = (char **)malloc(sizeof(char*)*MAX_LINE_PRINT);
		//entryid = (char **)malloc(sizeof(char*)*MAX_LINE_PRINT);
		//starts = (uint32_t *)malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
		starts = (uint32_t *)calloc(MAX_LINE_PRINT, sizeof(uint32_t));
		pends = (uint32_t *)malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
		values = (float *)malloc(sizeof(float) * MAX_LINE_PRINT);
		coverages = (uint16_t *)malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
		strands = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
		contexts = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
		int printL = 0, chrprinHdr = 0;
		fprintf(stderr, "[DM::calmeth] Processing chrom %d %d %d\n", genome_id, totalC, totalG);
		//--------DMC---------------//
		//
		int plus_mCG=0,plus_mCHG=0,plus_mCHH=0;
		int plusCG=0,plusCHG=0,plusCHH=0;
		int count_plus_CG=0,count_plus_CHG=0,count_plus_CHH=0;
		int Neg_mCG=0,Neg_mCHG=0,Neg_mCHH=0;
		int NegCG=0,NegCHG=0,NegCHH=0;
		int count_neg_CG=0,count_neg_CHG=0,count_neg_CHH=0;
		//------F
		
		int whichBins=0;int wBins=0;//which bins
		int i = genome_id;
		{
			//string Genome_Seq;Genome_Seq=args->Genome_List[i].Genome;

			// weight methylation level
		    int pluscountperCG=0,pluscountperCHG=0,pluscountperCHH=0,NegcountperCG=0,NegcountperCHG=0,NegcountperCHH=0,pluscountCG=0,NegcountCG=0,pluscountCHG=0,NegcountCHG=0,
		               pluscountCHH=0,NegcountCHH=0;
		    int pluscountperCG_1=0,pluscountperCHG_1=0,pluscountperCHH_1=0,NegcountperCG_1=0,NegcountperCHG_1=0,NegcountperCHH_1=0,pluscountCG_1=0,NegcountCG_1=0,pluscountCHG_1=0,
		            NegcountCHG_1=0,pluscountCHH_1=0,NegcountCHH_1=0;  
			int nb=0,nbins = ceil(args->Genome_Offsets[i].Offset/binspan);
			int nRegionBins = ceil(args->Genome_Offsets[i].Offset/RegionBins);
			char * Genome = args->Genome_Offsets[i].Genome;
			//fprintf(stderr, "chrom len %d\n", args->Genome_Offsets[i].Offset);
			for(int l=0;l<args->Genome_Offsets[i].Offset;l++)//loci
			{
				/*
		            if( (nb!=nbins && l==(nb+1)*binspan-1) || l==args->Genome_Offsets[i].Offset-1){
		                if(nb<=nbins && nb>0){
		                	if(pluscountCG+pluscountCG_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCG\n",args->Genome_Offsets[i].Genome,nb,((double)(pluscountperCG+pluscountperCG_1)/(pluscountCG+pluscountCG_1)) );
		                	if(pluscountCHG+pluscountCHG_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCHG\n",args->Genome_Offsets[i].Genome,nb,((double)(pluscountperCHG+pluscountperCHG_1)/(pluscountCHG+pluscountCHG_1)) );
		                	if(pluscountCHH+pluscountCHH_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCHH\n",args->Genome_Offsets[i].Genome,nb,((double)(pluscountperCHH+pluscountperCHH_1)/(pluscountCHH+pluscountCHH_1)) );

		                	if(NegcountCG+NegcountCG_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t-%f\tCG\n",args->Genome_Offsets[i].Genome,nb,((double)(NegcountperCG+NegcountperCG_1)/(NegcountCG+NegcountCG_1)));
		                	if(NegcountCHG+NegcountCHG_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t-%f\tCHG\n",args->Genome_Offsets[i].Genome,nb,((double)(NegcountperCHG+NegcountperCHG_1)/(NegcountCHG+NegcountCHG_1)));
		                	if(NegcountCHH+NegcountCHH_1>=Mcoverage) fprintf(BINsOUTFILE,"%s\t%d\t-%f\tCHH\n",args->Genome_Offsets[i].Genome,nb,((double)(NegcountperCHH+NegcountperCHH_1)/(NegcountCHH+NegcountCHH_1)) );
		                	
		                	if(pluscountCG+pluscountCG_1+NegcountCG+NegcountCG_1>coverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCG\n",Genome_Offsets[i].Genome,nb,((double)(pluscountperCG+pluscountperCG_1+NegcountperCG+NegcountperCG_1)/(pluscountCG+pluscountCG_1+NegcountCG+NegcountCG_1)) );
		                	if(pluscountCHG+pluscountCHG_1+NegcountCHG+NegcountCHG_1>coverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCHG\n",Genome_Offsets[i].Genome,nb,((double)(pluscountperCHG+pluscountperCHG_1+NegcountperCHG+NegcountperCHG_1)/(pluscountCHG+pluscountCHG_1+NegcountCHG+NegcountCHG_1)) );
		                	if(pluscountCHH+pluscountCHH_1+NegcountCHH+NegcountCHH_1>coverage) fprintf(BINsOUTFILE,"%s\t%d\t%f\tCHH\n",Genome_Offsets[i].Genome,nb,((double)(pluscountperCHH+pluscountperCHH_1+NegcountperCHH+NegcountperCHH_1)/(pluscountCHH+pluscountCHH_1+NegcountCHH+NegcountCHH_1)) );
		                } 
		                pluscountperCG_1=pluscountperCG;pluscountperCHG_1=pluscountperCHG;pluscountperCHH_1=pluscountperCHH;NegcountperCG_1=NegcountperCG;NegcountperCHG_1=NegcountperCHG;NegcountperCHH_1=NegcountperCHH;pluscountCG_1=pluscountCG;NegcountCG_1=NegcountCG;pluscountCHG_1=pluscountCHG;
		                NegcountCHG_1=NegcountCHG;pluscountCHH_1=pluscountCHH;NegcountCHH_1=NegcountCHH;
		                nb++;
		                pluscountperCG=0;pluscountperCHG=0;pluscountperCHH=0;NegcountperCG=0;NegcountperCHG=0;NegcountperCHH=0;pluscountCG=0;NegcountCG=0;pluscountCHG=0;NegcountCHG=0;pluscountCHH=0;NegcountCHH=0;
		            }
		            //-----------DMR region--//
				wBins=whichBins;
				whichBins = floor((double)l/RegionBins);
				if(wBins+1==whichBins || l==args->Genome_Offsets[i].Offset-1 )
				{//chr10   60441   -       CHH     0       4
					if(count_plus_CG>=nCs ) fprintf(REGION_OUT_CG,"%s\t%d\t+\tCG\t%d\t%d\n",Genome,wBins*RegionBins+1,plus_mCG,plusCG);
					if(count_neg_CG>=nCs) fprintf(REGION_OUT_CG,"%s\t%d\t-\tCG\t%d\t%d\n",Genome,wBins*RegionBins+1,Neg_mCG,NegCG);
					if(count_plus_CHG>=nCs ) fprintf(REGION_OUT_CHG,"%s\t%d\t+\tCHG\t%d\t%d\n",Genome,wBins*RegionBins+1,plus_mCHG,plusCHG);
					if(count_neg_CHG>=nCs) fprintf(REGION_OUT_CHG,"%s\t%d\t-\tCHG\t%d\t%d\n",Genome,wBins*RegionBins+1,Neg_mCHG,NegCHG);
					if(count_plus_CHH>=nCs ) fprintf(REGION_OUT_CHH,"%s\t%d\t+\tCHH\t%d\t%d\n",Genome,wBins*RegionBins+1,plus_mCHH,plusCHH);
					if(count_neg_CHH>=nCs) fprintf(REGION_OUT_CHH,"%s\t%d\t-\tCHH\t%d\t%d\n",Genome,wBins*RegionBins+1,Neg_mCHH,NegCHH);
					
					plus_mCG=0;plus_mCHG=0;plus_mCHH=0;plusCG=0;plusCHG=0;plusCHH=0;
					Neg_mCG=0;Neg_mCHG=0;Neg_mCHH=0;NegCG=0;NegCHG=0;NegCHH=0;
					count_plus_CG=0;count_plus_CHG=0;count_plus_CHH=0;
					count_neg_CG=0;count_neg_CHG=0;count_neg_CHH=0;
				}
				*/
				//-----F
				std::string context;
				
				if(args->Methy_List.plusMethylated[l]+args->Methy_List.plusUnMethylated[l]>=1) //coverage
				{//chromsome loci strand context methratio eff_CT_count C_count T_count CT_count rev_G_count rev_GA_count
					int C_count=args->Methy_List.plusMethylated[l];
					int T_count=args->Methy_List.plusUnMethylated[l];
					int rev_G=args->Methy_List.NegG[l];
					int rev_A=args->Methy_List.NegA[l];
					float revGA=0;
					if( (rev_G+rev_A)> 20 ) revGA=float(rev_G)/float(rev_G+rev_A); 
					//context
					
					string Fivecontext;
					//char genome_Char = toupper(args->Genome_List[i].Genome[l]);
					if(l+2+1 < args->Genome_Offsets[i].Offset)
					{
						if(l>=2 ) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 5); //Genome_Seq.substr(l-2,5);
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 4); //Genome_Seq.substr(l-1,4);
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 3); //Genome_Seq.substr(l,3);
					}else if(l+1+1 < args->Genome_Offsets[i].Offset)
					{
						if(l>=2 ) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 4)+"N";
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 3)+"N";
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 2)+"N";
					}else
					{
						if(l>=2 ) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 3)+"NN";
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 2)+"NN";
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 1)+"NN";
					}
					int stringlength=Fivecontext.length();
											
					transform(Fivecontext.begin(), Fivecontext.end(), Fivecontext.begin(), ::toupper);

					char charFor1='N',charFor2='N';//printf("\n%s\n",Fivecontext.c_str());
					if(l+1< args->Genome_Offsets[i].Offset) charFor1=toupper(args->Genome_List[i].Genome[l+1]);
					if(l+2< args->Genome_Offsets[i].Offset) charFor2=toupper(args->Genome_List[i].Genome[l+2]);
					int Cover = C_count+T_count;
					for(int nc=0; nc<15; nc++){
						if(Cover>nc){ CoverC[nc]++; }
					}

					if(charFor1=='G') //(args.Methy_List[i].MethContext[l]==1)
					{
						context="CG";
						plus_mCGcount+=C_count;
						plusCGcount+=(C_count+T_count);
						//--DMR
						plus_mCG+=C_count;
						plusCG+=(C_count+T_count);
						count_plus_CG++;
						//chromsome bins
						pluscountperCG+=C_count;
						pluscountCG+=(C_count+T_count);
					}
					else if(charFor1!='G' && charFor2=='G') //(args.Methy_List[i].MethContext[l]==2)
					{
						context="CHG";
						plus_mCHGcount+=C_count;
						plusCHGcount+=(C_count+T_count);
						//--DMR
						plus_mCHG+=C_count;
						plusCHG+=(C_count+T_count);
						count_plus_CHG++;
						//bins
						pluscountperCHG+=C_count;
						pluscountCHG+=(C_count+T_count);								
					}
					else if(charFor1!='G' && charFor1!='G') //(args.Methy_List[i].MethContext[l]==3)
					{
						context="CHH";
						plus_mCHHcount+=C_count;
						plusCHHcount+=(C_count+T_count);
						//DMR
						plus_mCHH+=C_count;
						plusCHH+=(C_count+T_count);
						count_plus_CHH++;
						//bins
						pluscountperCHH+=C_count;
						pluscountCHH+=(C_count+T_count);							
					}
					else context="NA";
					
					if(Cover<Mcoverage || Cover > maxcoverage)
						continue;
					//methratio	
					float PlusMethratio;
					if(revGA>0) 
						PlusMethratio=std::min(float(C_count)/(float(C_count+T_count)*revGA),(float)1.0);
					else 
						PlusMethratio=float(C_count)/float(C_count+T_count);
					
					unsigned int mCdensityloci=std::min( (int)floor((double)(PlusMethratio*100) )  ,99);
					if( !strcmp(context.c_str(),"CG")) mCGdensity[mCdensityloci]++;
					else if(!strcmp(context.c_str(),"CHG")) mCHGdensity[mCdensityloci]++;
					else if(!strcmp(context.c_str(),"CHH")) mCHHdensity[mCdensityloci]++;							
					
					string category = "NA";
						if(PlusMethratio>=0.8){
							M++;
							category="M";
							if(!strcmp(context.c_str(),"CG")) M_CG++;
						}else if(PlusMethratio >=0.6){
							Mh++;
							category="Mh";
							if(!strcmp(context.c_str(),"CG")) Mh_CG++;
						}else if(PlusMethratio >=0.4){
							H_AllC++;
							category="H";
							if(!strcmp(context.c_str(),"CG")) H_CG++;
						}else if(PlusMethratio >=0.2){
							hU++;
							category="hU";
							if(!strcmp(context.c_str(),"CG")) hU_CG++;
						}else{
							U++;
							category="U";
							if(!strcmp(context.c_str(),"CG")) U_CG++;
						}

					//print store
					chromsUse[printL] = strdup(args->Genome_Offsets[i].Genome);
					starts[printL] = l+1;
					pends[printL] = l+2;
					coverages[printL] = (C_count+T_count);
                	values[printL] = PlusMethratio;
					strands[printL] = 0; //0 represent '+'
					if(strcmp(context.c_str(), "C") == 0){
						contexts[printL] = 0;
					}else if(strcmp(context.c_str(), "CG") == 0){
						contexts[printL] = 1;
					}else if(strcmp(context.c_str(), "CHG") == 0){
						contexts[printL] = 2;
					}else{ // if(strcmp(context.c_str(), "CHH") == 0){
						contexts[printL] = 3;
					}
					printL++;
					if(chrprinHdr==0){
						int response = bmAddIntervals(fp, chromsUse, starts, pends, values, coverages, strands, contexts, 
		                entryid, 1);
        		        if(response) goto error;
					}
					else if(printL>MAX_LINE_PRINT){
						//We can continue appending similarly formatted entries
						//N.B. you can't append a different chromosome (those always go into different
						if(bmAppendIntervals(fp, starts+1, pends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) goto error;
						printL = 0;
					}
					chrprinHdr = 1;

					if(printtxt == 1){  
						if(revGA>0) fprintf(METHOUTFILE,"%s\t%d\t+\t%s\t%d\t%d\t%f\t%0.001f\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),PlusMethratio,float(C_count+T_count)*revGA,rev_G,(rev_A+rev_G),category.c_str(),Fivecontext.c_str());
						else fprintf(METHOUTFILE,"%s\t%d\t+\t%s\t%d\t%d\t%f\tnull\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),PlusMethratio,rev_G,(rev_A+rev_G),category.c_str(),Fivecontext.c_str());
					}

					/*
					if(!strcmp(context.c_str(),"CG")) fprintf(LOC_OUT_CG,"%s\t%d\t+\tCG\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					else if(!strcmp(context.c_str(),"CHG")) fprintf(LOC_OUT_CHG,"%s\t%d\t+\tCHG\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					else if(!strcmp(context.c_str(),"CHH")) fprintf(LOC_OUT_CHH,"%s\t%d\t+\tCHH\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					*/
				}

				if(args->Methy_List.NegMethylated[l]+args->Methy_List.NegUnMethylated[l]>=1) //coverage
				{
					int C_count=args->Methy_List.NegMethylated[l];
					int T_count=args->Methy_List.NegUnMethylated[l];
					int rev_G=args->Methy_List.plusG[l];
					int rev_A=args->Methy_List.plusA[l];
					float revGA=0;
					if( (rev_G+rev_A)>20 ) revGA=float(rev_G)/float(rev_G+rev_A); 
					//context
					string Fivecontext;
					//char genome_Char = toupper(Genome_Seq[l]);
					if(l+2+1 < args->Genome_Offsets[i].Offset)
					{
						if(l>=2 ) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 5);
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 4);
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 3);
					}else if(l+1+1 < args->Genome_Offsets[i].Offset)
					{
						if(l>=2) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 4)+"N";
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 3)+"N";
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 2)+"N";
					}else{
						if(l>=2 ) Fivecontext= getstring(args->Genome_List[i].Genome, l-2, 3)+"NN";
						else if(l==1) Fivecontext = "N" + getstring(args->Genome_List[i].Genome, l-1, 2)+"NN";
						else if(l==0) Fivecontext = "NN" + getstring(args->Genome_List[i].Genome, l, 1)+"NN";
					}
					
					int stringlength=strlen(Fivecontext.c_str());
					char Fcontext[6];
					memcpy(Fcontext,Fivecontext.c_str(),stringlength+1);
					if(stringlength<5) printf("\n%d %d %d %s \n",strlen(Fivecontext.c_str()), l, args->Genome_Offsets[i].Offset, args->Genome_Offsets[i].Genome);
					ReverseC_Context(Fcontext,Fivecontext.c_str(),stringlength);
					
					int Cover = C_count+T_count;
					for(int nc=0; nc<15; nc++){
						if(Cover>nc){ CoverC[nc]++; }
					}

					char charBac1='N',charBac2='N';
					if(l>=1) charBac1=toupper(args->Genome_List[i].Genome[l-1]);
					if(l>=2) charBac2=toupper(args->Genome_List[i].Genome[l-2]);
					if(charBac1=='C') //(args.Methy_List[i].MethContext[l]==1)
					{
						context="CG";
						Neg_mCGcount+=C_count;
						NegCGcount+=(C_count+T_count);
						//--DMR
						Neg_mCG+=C_count;
						NegCG+=(C_count+T_count);
						count_neg_CG++;
						//bins
		                          NegcountperCG+=C_count;
		                          NegcountCG+=(C_count+T_count);
					}
					else if(charBac1!='C' && charBac2=='C') //(args.Methy_List[i].MethContext[l]==2)
					{
						context="CHG";
						Neg_mCHGcount+=C_count;
						NegCHGcount+=(C_count+T_count);
						//---DMR
						Neg_mCHG+=C_count;
						NegCHG+=(C_count+T_count);
						count_neg_CHG++;
						//bins
		                        NegcountperCHG+=C_count;
		                        NegcountCHG+=(C_count+T_count);
					}
					else if(charBac1!='C' && charBac2!='C') //(args.Methy_List[i].MethContext[l]==3)
					{
						context="CHH";
						Neg_mCHHcount+=C_count;
						NegCHHcount+=(C_count+T_count);
						//---DMR
						Neg_mCHH+=C_count;
						NegCHH+=(C_count+T_count);
						count_neg_CHH++;
						//bins
		                        NegcountperCHH+=C_count;
		                        NegcountCHH+=(C_count+T_count);
					}
					else context="NA";

                    if(Cover<Mcoverage || Cover > maxcoverage) continue;
					//methratio	
					float NegMethratio;
					if(revGA>0)
						NegMethratio=std::min(float(C_count)/(float(C_count+T_count)*revGA),(float)1.0);
					else 
						NegMethratio=float(C_count)/(float(C_count+T_count));
					
					unsigned int mCdensityloci= std::min( (int)floor((double)(NegMethratio*100) ) ,99);
					if(!strcmp(context.c_str(),"CG")) mCGdensity[mCdensityloci]++;
					else if(!strcmp(context.c_str(),"CHG")) mCHGdensity[mCdensityloci]++;
					else if(!strcmp(context.c_str(),"CHH")) mCHHdensity[mCdensityloci]++;
					
					string category = "NA";
			              if(NegMethratio>=0.8){
			                  M++;
			                  category="M";
			                  if(!strcmp(context.c_str(),"CG")) M_CG++;
			              }else if(NegMethratio >=0.6){
			                  Mh++;
			                  category="Mh";
			                  if(!strcmp(context.c_str(),"CG")) Mh_CG++;
			              }else if(NegMethratio >=0.4){
			                  H_AllC++;
			                  category="H";
			                  if(!strcmp(context.c_str(),"CG")) H_CG++;
			              }else if(NegMethratio >=0.2){
			                  hU++;
			                  category="hU";
			                  if(!strcmp(context.c_str(),"CG")) hU_CG++;
			              }else{
			                  U++;
			                  category="U";
			                  if(!strcmp(context.c_str(),"CG")) U_CG++;
			              }
					//print store
					chromsUse[printL] = strdup(args->Genome_Offsets[i].Genome);
					starts[printL] = l+1;
					pends[printL] = l+2;
					coverages[printL] = (C_count+T_count);
                	values[printL] = NegMethratio;
					strands[printL] = 1; //1 represent '-'
					if(strcmp(context.c_str(), "C") == 0){
						contexts[printL] = 0;
					}else if(strcmp(context.c_str(), "CG") == 0){
						contexts[printL] = 1;
					}else if(strcmp(context.c_str(), "CHG") == 0){
						contexts[printL] = 2;
					}else{ // if(strcmp(context.c_str(), "CHH") == 0){
						contexts[printL] = 3;
					}
					printL++;
					if(chrprinHdr==0){
						int response = bmAddIntervals(fp, chromsUse, starts, pends, values, coverages, strands, contexts, 
		                entryid, 1);
        		        if(response) goto error;
					}
					else if(printL>MAX_LINE_PRINT){
						//We can continue appending similarly formatted entries
						//N.B. you can't append a different chromosome (those always go into different
						if(bmAppendIntervals(fp, starts+1, pends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) goto error;
						printL = 0;
					}
					chrprinHdr = 1;

					if(printtxt == 1){
						if(revGA>0) fprintf(METHOUTFILE,"%s\t%d\t-\t%s\t%d\t%d\t%f\t%0.001f\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),NegMethratio,float(C_count+T_count)*revGA,rev_G,(rev_G+rev_A),category.c_str(),Fcontext);
						else fprintf(METHOUTFILE,"%s\t%d\t-\t%s\t%d\t%d\t%f\tnull\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),NegMethratio,rev_G,(rev_G+rev_A),category.c_str(),Fcontext);
					}
					/*
					if(!strcmp(context.c_str(),"CG")) fprintf(LOC_OUT_CG,"%s\t%d\t-\tCG\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					else if(!strcmp(context.c_str(),"CHG")) fprintf(LOC_OUT_CHG,"%s\t%d\t-\tCHG\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					else if(!strcmp(context.c_str(),"CHH")) fprintf(LOC_OUT_CHH,"%s\t%d\t-\tCHH\t%d\t%d\n",Genome,l+1,C_count,(C_count+T_count));
					*/
				}
			}
		}
		//end print
		if(printL > 1) {
            fprintf(stderr, "[DM::calmeth] print %d %d %d\n", starts[printL-1], pends[printL-1], printL-1);
            if(bmAppendIntervals(fp, starts+1, pends+1, values+1, coverages+1, strands+1, contexts+1, entryid, printL-1)) goto error;
            printL = 0; 
        }
		//
		fprintf(stderr, "[DM::calmeth] Free mem in Methy_List\n");
		for(i=0; i< longestChr; i++){
			args->Methy_List.plusG[i] = 0;
			args->Methy_List.plusA[i] = 0;
			args->Methy_List.NegG[i] = 0;
			args->Methy_List.NegA[i] = 0;
			args->Methy_List.plusMethylated[i] = 0;
			args->Methy_List.plusUnMethylated[i] = 0;
			args->Methy_List.NegMethylated[i] = 0;
			args->Methy_List.NegUnMethylated[i] = 0;
		}
		//printf("\n");

		fprintf(stderr, "[DM::calmeth] Free mem in chromsUse\n");
		for(i =0; i < MAX_LINE_PRINT; i++){
            if(starts[i]>0 && chromsUse[i]) free(chromsUse[i]);
			//fprintf(stderr, "Free mem in chromsUse %d\n", i);
			// //if(entryid[i]) free(entryid[i]);
			//fprintf(stderr, "Free mem in chromsUse2 %d\n", i);
        }
		fprintf(stderr, "[DM::calmeth] Free mem in all others\n\n");
        free(chromsUse); //free(entryid);

        free(starts);
        free(pends); free(values); free(coverages); free(strands); free(contexts);

	}//end methratio
	return;

error:
	fprintf(stderr, "\nEEEEE main Received an error somewhere!\n");
}

char *fastStrcat(char *s, char *t)
{
	assert(s != NULL && t != NULL);

	while (*s != '\0')
		s++;

	//while ((*s++ = *t++) != '\0');
	while(1){
		if(*t == '\0' || *t == '\n' || *t == '\r'){
			*s++ = '\0';
			break;
		}
		*s++ = *t++;
	}

	return --s;
}

int conutMismatch(string CIGr, int chrLen, char* Genome_seq, string readString, int pos, int hitType)
{
	char temp[5];unsigned lens=0;int Glens=0;int RLens=0;
	unsigned n=0;
	int Nmismatch=0;  //chrLen; //((ARGS *)arg)->Genome_Offsets[Hash_Index+1].Offset
	const char* cigr=CIGr.c_str();
	while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
	{
		if(*cigr>='0' && *cigr<='9')
		{
			temp[n]=*cigr;
			cigr++;n++;
		}else if(*cigr=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			lens+=length;
	        	RLens+=length;
        	    	cigr++;n=0;
		}else if(*cigr=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			for(int k=lens,r=RLens,g=Glens;k<length+lens;r++,k++,g++)
			{
				if (pos+g-1 >= chrLen) break;

				char genome_Char = toupper(Genome_seq[pos+g-1]);
				if (hitType==1 || hitType==3) {
					if (readString[k] != genome_Char && !(readString[k]=='T' && genome_Char=='C')) 
					{
						Nmismatch++;
					}
				}
				else if (hitType==2 || hitType==4) {
					if (readString[k] != genome_Char && !(readString[k]=='A' && genome_Char=='G')) {
						
						Nmismatch++;
					}
				}
			}
			cigr++;n=0;lens+=length;Glens+=length;RLens+=length;
		}else if(*cigr=='I')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			lens+=length;
			RLens+=length;
			cigr++;n=0;
		}else if(*cigr=='D')
		{
			int i;temp[n]='\0';unsigned int length=atoi(temp);
			Glens+=length;RLens+=length;
			cigr++;n=0;
		}else
		{
			break;
		}
	}
	return Nmismatch;
}

int processbamread(const bam_hdr_t *header, const bam1_t *b, char* Dummy,int &Flag,char* Chrom,int &pos,int &mapQuality,char* CIG,char* Chrom_P,int &pos_P,int &Insert_Size,char* forReadString,char* forQuality, int &hitType)
{
        uint8_t *s = bam_get_seq(b), *t = bam_get_qual(b);
        int i;
        const bam1_core_t *c = &b->core;
	strcpy(Dummy,  bam_get_qname(b));
	Flag = c->flag;

	if( (Flag & 0x100) || (Flag & 0x200) || (Flag & 0x400) || (Flag & 0x800) || (Flag & 0x4))
        	return -1;
	if( !(Flag & 0x1) )
        {
        	if(Flag==0)
                	hitType=1;
                else if(Flag==16)
                	hitType=4;
        }else if( !(Flag & 0x10)  && (Flag & 0x40) )
        	hitType=1;
       	else if( !(Flag & 0x10)  && (Flag & 0x80) )
        	hitType=2;
        else if( (Flag & 0x10)  && (Flag & 0x80) )
        	hitType=3;
        else if( (Flag & 0x10)  && (Flag & 0x40) )
        	hitType=4;
        if(hitType==0) return -1;


        if (c->tid < 0) return -1;
        else {
                if (header) strcpy(Chrom, header->target_name[c->tid]); 
                else sprintf(Chrom, "%d", c->tid); 
        }
	pos = c->pos + 1;
	mapQuality = c->qual;

	//define
	if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) return -1;
	char strtemp[256];int j=0;
        if (c->n_cigar == 0) return -1;
        else {
                for (i = 0; i < c->n_cigar; ++i) {
                    	sprintf(strtemp, "%d%c", bam_get_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP=X"[bam_get_cigar(b)[i]&BAM_CIGAR_MASK]);
			for(int l= 0; l<strlen(strtemp); ++l,j++) CIG[j] = strtemp[l];
                }
        }
	CIG[j] = '\0';
        if (c->mtid < 0) Chrom_P[0] =  '*';
        else if (c->mtid == c->tid) sprintf(Chrom_P, "="); 
        else {
                if (header) strcpy(Chrom_P, header->target_name[c->mtid]); 
                else sprintf(Chrom_P, "%d", c->mtid); 
        }
	if(strcmp(Chrom_P, "*") == 0) {
        	pos_P = 0;
                Insert_Size = 0;
        }else{
		pos_P = c->mpos + 1;
		Insert_Size=c->isize;
	}
        if (c->l_qseq) {
            for (i = 0; i < c->l_qseq; ++i) forReadString[i] = seq_nt16_str[bam_seqi(s, i)];
		    forReadString[i] = '\0';
            if (t[0] == 0xff) forQuality[0] =  '*';
            else for (i = 0; i < c->l_qseq; ++i) forQuality[i] = (char)(t[i] + 33); 
		    if(i!=0) forQuality[i] = '\0';
        } else {forReadString[0] = '*'; forQuality[0] =  '*';}
	
	return 1;

        s = bam_get_aux(b);
	char read[100];
        while (s < b->data + b->l_data) {
                uint8_t type, key[2];
                key[0] = s[0]; key[1] = s[1];
                s += 2; type = *s; ++s;
                sprintf(read, "\t%s:", (char*)key); 
                if (type == 'A') { sprintf(read, "A:%c", *s); ++s; }
                else if (type == 'C') { sprintf(read, "i:%d", *s); ++s; }
                else if (type == 'c') { sprintf(read, "i:%d", *(int8_t*)s); ++s; }
                else if (type == 'S') { sprintf(read, "i:%d", *(uint16_t*)s); s += 2; }
                else if (type == 's') { sprintf(read, "i:%d", *(int16_t*)s); s += 2; }
                else if (type == 'I') { sprintf(read, "i:%d", *(uint32_t*)s); s += 4; }
                else if (type == 'i') { sprintf(read, "i:%d", *(int32_t*)s); s += 4; }
                else if (type == 'f') { sprintf(read, "f:%g", *(float*)s); s += 4; }
                else if (type == 'd') { sprintf(read, "d:%lg", *(double*)s); s += 8; }
                else if (type == 'Z' || type == 'H') { sprintf(read, "%c:", type); while (*s) sprintf(read, "%c", *s++); ++s; }
                else if (type == 'B') {
                        uint8_t sub_type = *(s++);
                        int32_t n;
                        memcpy(&n, s, 4);
                        s += 4; // no point to the start of the array
                    	sprintf(read, "%c:%c", type, sub_type);
                        for (i = 0; i < n; ++i) {
                                sprintf(read,",");
                                if ('c' == sub_type || 'c' == sub_type) { sprintf(read, "%d", *(int8_t*)s); ++s; }
                                else if ('C' == sub_type) { sprintf(read, "%d", *(uint8_t*)s); ++s; }
                                else if ('s' == sub_type) { sprintf(read, "%d", *(int16_t*)s); s += 2; }
                                else if ('S' == sub_type) { sprintf(read, "%d", *(uint16_t*)s); s += 2; }
                                else if ('i' == sub_type) { sprintf(read, "%d", *(int32_t*)s); s += 4; }
                                else if ('I' == sub_type) { sprintf(read, "%d", *(uint32_t*)s); s += 4; }
                                else if ('f' == sub_type) { sprintf(read, "%g", *(float*)s); s += 4; }
                        }
                }
        }
        
}

void *Process_read(void *arg)
{
	//unsigned Total_Reads=0, Total_mapped = 0, forward_mapped = 0, reverse_mapped = 0;
	string fileIS = Prefix + ".insert_size.1M.txt";
	long process_vali = 0;
	FILE* fIS = File_Open(fileIS.c_str(), "w");
	int processed_vali = 0;
	int Progress=0;Number_of_Tags=INITIAL_PROGRESS_READS;
	Init_Progress();
	int mismatch=0;int pos=0;int Top_Penalty=0;int mapQuality=0;int Flag=-1;
	string readString="";
	int hitType=0;
    int fileprocess = 0;
	
	string hits[MAX_HITS_ALLOWED];
	char Comp_String[MAXTAG];for (int i=1;i<MAXTAG;Comp_String[i++]=0);
	//start to read batman hit file
	char *s2t = (char*) malloc(1000);
	char read_Methyl_Info[600];char rawReadBeforeBS[600];char temp[5];
	char Dummy[BATBUF],forReadString[BATBUF],Chrom[CHROMSIZE];
	char Chrom_P[CHROMSIZE];int pos_P=0;int Insert_Size=0;int Qsingle=0; //Paired-end reads
	string CIGr;char CIG[BATBUF];
	char forQuality[BATBUF],rcQuality[BATBUF],Quality[BATBUF];
	int r=1;char* samaddress;
    struct timeval now;
    struct timespec outtime;
	bam1_t *b = bam_init1();
    hts_idx_t *idx = NULL;
    hts_itr_t *iter;
    if(strcmp(((ARGS *)arg)->processChr, "NAN-mm") != 0){
        if ((idx = sam_index_load(((ARGS *)arg)->BamInFile, ((ARGS *)arg)->INbamfilename )) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            exit(0);
        }
        if ((iter = sam_itr_querys(idx, ((ARGS *)arg)->header, ((ARGS *)arg)->processChr)) == 0) {
            fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, ((ARGS *)arg)->processChr);
            exit(0);
        }
    }
	while( (!bamformat && (samaddress = fgets(s2t,BATBUF,((ARGS *)arg)->samINFILE))!=NULL) || (bamformat && r>0 ))
	{
		if(r < -1) {
			fprintf(stderr, "\ntruncated file.\n");
		}
		hitType = 0;
		
		Progress++;
		fileprocess++;

		if(bamformat) 
		{
            if(strcmp(((ARGS *)arg)->processChr, "NAN-mm") == 0){
			    //(r = samread(( (ARGS *)arg)->BamInFile, b));
                r = sam_read1(( (ARGS *)arg)->BamInFile, ((ARGS *)arg)->header, b);
			    //bam_tostring(((ARGS *)arg)->header , b, s2t);
			    int ct = processbamread(((ARGS *)arg)->header, b, Dummy,Flag,Chrom,pos,mapQuality,CIG,Chrom_P,pos_P,Insert_Size,forReadString,forQuality, hitType);
			    if(ct == -1) continue;
			    if(r<=0) break;
            }else{
                r = sam_itr_next(((ARGS *)arg)->BamInFile, iter, b);
                int ct = processbamread(((ARGS *)arg)->header, b, Dummy,Flag,Chrom,pos,mapQuality,CIG,Chrom_P,pos_P,Insert_Size,forReadString,forQuality, hitType);
                if(ct == -1) continue;
                if(r<=0) break;
            }
		}
        Total_Reads++;
		
		if ( fileprocess>=1000000  ) {
			fprintf_time(stderr, "Processed %d reads.\n", Total_Reads);
			fileprocess = 0;
        }

		if(s2t[0]=='@') 
		{
			continue;
		}

        printheader = false;
		if(!bamformat)
			sscanf(s2t,"%s%d%s%d%d%s%s%d%d%s%s",Dummy,&Flag,Chrom,&pos,&mapQuality,CIG,Chrom_P,&pos_P,&Insert_Size,forReadString,forQuality);

		map<string, int>::iterator iter;
		int H = -1;
		if(strcmp(newchr.c_str(), Chrom) != 0 ) newchr = Chrom;
		
		if(strcmp(processingchr.c_str(), "NULL") == 0 ) processingchr = Chrom;
		if(strcmp(processingchr.c_str(), Chrom) != 0) {
            fprintf(stderr, "[DM::calmeth] Print output of %s\n", processingchr.c_str());
            iter = String_Hash.find(processingchr.c_str());
			if(iter != String_Hash.end()){
				H = iter->second;
				print_meth_tofile(H, (ARGS *)arg);
			}
	        processingchr = Chrom;
		}

		if(!bamformat){
 			if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) continue;
			if(strcmp(Chrom_P, "*") == 0) {
				pos_P = 0;
				Insert_Size = 0;
			}
		}

		readString=forReadString;
		int Read_Len=readString.length();
		CIGr=CIG;
		//for(;forReadString[Read_Len]!=0 && forReadString[Read_Len]!='\n' && forReadString[Read_Len]!='\r';Read_Len++);
	    	iter = String_Hash.find(processingchr.c_str());
		H = -1;
		if(iter != String_Hash.end()){
			H = iter->second;
		}else continue;

		if(Flag==-1) printf("\n%s\n", Dummy);
		assert(Flag!=-1);
		if(!bamformat){
			if( (Flag & 0x100) || (Flag & 0x200) || (Flag & 0x400) || (Flag & 0x800) || (Flag & 0x4))
        	    continue;
			if( !(Flag & 0x1) )
			{
				if(Flag==0)
					hitType=1;
				else if(Flag==16)
					hitType=4;
			}else if( !(Flag & 0x10)  && (Flag & 0x40) )
				hitType=1;
			else if( !(Flag & 0x10)  && (Flag & 0x80) )
				hitType=2;
			else if( (Flag & 0x10)  && (Flag & 0x80) )
				hitType=3;
			else if( (Flag & 0x10)  && (Flag & 0x40) )
				hitType=4;
			if(hitType==0) continue;
		}
		// !(Flag & 0x100) && !(Flag & 0x200)
		Total_mapped++;
		if(hitType == 1 || hitType == 2) forward_mapped++;
		else if(hitType == 4 || hitType == 3) reverse_mapped++;
		if(strcmp(Chrom_P, "=") == 0) {
			if(process_vali <= 1000000){
				if(Insert_Size>=0 && Insert_Size<=2000){
					fprintf(fIS, "%d\n", Insert_Size);
					process_vali++;
				}
			}
		}

		int Nmismatch=conutMismatch(CIGr, ((ARGS *)arg)->Genome_Offsets[H].Offset, ((ARGS *)arg)->Genome_List[H].Genome, readString, pos, hitType);
		if(Nmismatch > 0.5 + UPPER_MAX_MISMATCH * strlen(readString.c_str())) continue;
		//char* Genome;
		//
		//Genome=((ARGS *)arg)->Genome_List[H].Genome;//load current genome..
		
		unsigned long G_Skip=((ARGS *)arg)->Genome_List[H].Genome - ((ARGS *)arg)->Org_Genome;
            int Flag_rm=0;
		if(hitType == 1 || hitType == 3 ) Flag_rm=2; else Flag_rm=4;
		char Mark = '0';char MarkE='0';
		if(REMOVE_DUP){
            Mark=((ARGS *)arg)->Marked_Genome[pos+G_Skip];
            MarkE=((ARGS *)arg)->Marked_GenomeE[pos+G_Skip+readString.size()];
        }
        if( !REMOVE_DUP || (!Mark || !(Mark & Flag_rm)) || (!MarkE || !(MarkE & Flag_rm)) )
		{
			int Hash_Index=((ARGS *)arg)->Genome_List[H].Index;//load current genome..
			strcpy(rawReadBeforeBS,readString.c_str());
			unsigned lens=0;int Glens=0;int RLens=0;
			unsigned n=0;bool CONTINUE=false;
			const char* cigr=CIGr.c_str();
			while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
			{
				if(*cigr>='0' && *cigr<='9')
				{
					temp[n]=*cigr;
					cigr++;n++;
				}else if(*cigr=='S')
				{
					int i;temp[n]='\0';int length=atoi(temp);
					for(i=RLens;i<RLens+length;i++)
					{
						read_Methyl_Info[i] = 'S';
					}
					lens+=length;
	                RLens+=length;
	                cigr++;n=0;
				}else if(*cigr=='M')
				{
					temp[n]='\0';int length=atoi(temp);
					for(int k=lens,r=RLens,g=Glens;k<length+lens;r++,k++,g++)
					{
						read_Methyl_Info[r] = '=';
						if (pos+g-1 >= ((ARGS *)arg)->Genome_Offsets[Hash_Index].Offset) break;

						char genome_Char = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-1]);//
						char genome_CharFor1 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g+1-1]);
						char genome_CharFor2 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g+2-1]);
						char genome_CharBac1,genome_CharBac2;
						if(pos+g-1 > 2)
						{
							genome_CharBac1 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-1-1]);
							genome_CharBac2 = toupper(((ARGS *)arg)->Genome_List[H].Genome[pos+g-2-1]);
						}
						if (hitType==1 || hitType==3) {
							//pthread_mutex_lock(&meth_counter_mutex13);
							if (readString[k]=='C' && genome_Char=='C')
							{
								read_Methyl_Info[r] = 'M';
								if(Methratio ) ((ARGS *)arg)->Methy_List.plusMethylated[pos+g-1]++;
//if(pos+g-1 < 6) printf("\n%d %c\n", pos+g-1, Genome[pos+g-1]);
							//	if(hitType==1) {
									if(genome_CharFor1=='G') 
									{
										read_Methyl_Info[r] = 'Z';met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}//Z methylated C in CpG context
									else if(genome_CharFor1!='G' && genome_CharFor2!='G')
									{
										read_Methyl_Info[r] = 'H';met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}//H methylated C in CHH context
									else if(genome_CharFor1!='G' && genome_CharFor2=='G')
									{
										read_Methyl_Info[r] = 'X';met_CHG++;
										//if(Methratio && hitsToken[ind][4] > QualCut  )  ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}//X methylated C in CHG context
							}
							else if (readString[k]=='T' && genome_Char=='C')
							{
								read_Methyl_Info[r] = 'U';
								rawReadBeforeBS[k] = 'C';
								if(Methratio ) ((ARGS *)arg)->Methy_List.plusUnMethylated[pos+g-1]++;
							//	if(hitType==1) {
									if(genome_CharFor1=='G')
									{
										read_Methyl_Info[r] = 'z';non_met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}//z unmethylated C in CpG context
									else if(genome_CharFor1!='G' && genome_CharFor2!='G') 
									{
										read_Methyl_Info[r] = 'h';non_met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}//h unmethylated C in CHH context
									else if(genome_CharFor1!='G' && genome_CharFor2=='G') 
									{
										read_Methyl_Info[r] = 'x';non_met_CHG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}//x unmethylated C in CHG context
								
							}
							else if (readString[k] != genome_Char) 
							{
								read_Methyl_Info[r] = readString[k]; // genome_Char; for hypol
							}
							if(Methratio)
							{
								if(readString[k]=='G') ((ARGS *)arg)->Methy_List.plusG[pos+g-1]++;
								else if(readString[k]=='A') ((ARGS *)arg)->Methy_List.plusA[pos+g-1]++;
								//Methy_List[H].plusCover[pos+g-1]++;
							}
							//pthread_mutex_unlock(&meth_counter_mutex13);
						}
						else if (hitType==2 || hitType==4) {
							//pthread_mutex_lock(&meth_counter_mutex24);
							if (readString[k]=='G' && genome_Char=='G')
							{
								read_Methyl_Info[r] = 'M';
								if(Methratio ) ((ARGS *)arg)->Methy_List.NegMethylated[pos+g-1]++;
							//	if(hitType==2) 
							//	{
									if(genome_CharBac1=='C') 
									{
										read_Methyl_Info[r] = 'Z';met_CG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') 
									{
										read_Methyl_Info[r] = 'H';met_CHH++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') 
									{
										read_Methyl_Info[r] = 'X';met_CHG++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}
							}
							else if (readString[k]=='A' && genome_Char=='G')
							{
								read_Methyl_Info[r] = 'U';
								rawReadBeforeBS[k] = 'G';
								if(Methratio) ((ARGS *)arg)->Methy_List.NegUnMethylated[pos+g-1]++;
							//	if(hitType==2) 
							//	{
									if(genome_CharBac1=='C') 
									{
										read_Methyl_Info[r] = 'z';non_met_CG++;
										//if(Methratio && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=1;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') 
									{
										read_Methyl_Info[r] = 'h';non_met_CHH++;
										//if(Methratio && hitsToken[ind][4] > QualCut ) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=3;
									}
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') 
									{
										read_Methyl_Info[r] = 'x';non_met_CHG++;
										//if(Methratio  && hitsToken[ind][4] > QualCut) ((ARGS *)arg)->Methy_List[H].MethContext[pos+g-1]=2;
									}
							}
							else if (readString[k] != genome_Char) {
								
								read_Methyl_Info[r] = readString[k]; //genome_Char;
							}
							if(Methratio)
							{
								if(readString[k]=='C') ((ARGS *)arg)->Methy_List.NegG[pos+g-1]++;
								else if(readString[k]=='T') ((ARGS *)arg)->Methy_List.NegA[pos+g-1]++;
								//Methy_List[H].NegCover[pos+g-1]++;
							}
							//pthread_mutex_unlock(&meth_counter_mutex24);
						}
					}
					cigr++;n=0;lens+=length;Glens+=length;RLens+=length;
				}else if(*cigr=='I')
				{
					int i;temp[n]='\0';int length=atoi(temp);
					for(i=0;i<length;i++)
					{
						read_Methyl_Info[i+RLens] = 'I';
					}
					lens+=length;
					RLens+=length;
					cigr++;n=0;
				}else if(*cigr=='D')
				{
					int i;temp[n]='\0';unsigned int length=atoi(temp);
					for(i=0;i<length;i++)
					{
						read_Methyl_Info[i+RLens] = 'D';
					}
					Glens+=length;RLens+=length;
					cigr++;n=0;
				}else
				{
					//printf(" --%d%c-- ",atoi(temp),*cigr);
					CONTINUE=true;
					break;
					//continue;
					//exit(0);
				}
			}
			if(CONTINUE) continue;
							
			read_Methyl_Info[RLens]='\0';
			rawReadBeforeBS[lens]='\0';
			string mappingStrand="N";
			if(hitType==1) mappingStrand="YC:Z:CT\tYD:Z:f";
			else if(hitType==4) mappingStrand="YC:Z:CT\tYD:Z:r";
			else if(hitType==3) mappingStrand="YC:Z:GA\tYD:Z:r";
			else if(hitType==2) mappingStrand="YC:Z:GA\tYD:Z:f";
            if(Sam && ((ARGS *)arg)->OUTFILE != NULL) {
			    flockfile(((ARGS *)arg)->OUTFILE);
			    if(Sam ) //&& Nmismatch <= UPPER_MAX_MISMATCH )
			    {
				if(SamSeqBeforeBS) 
				{
					fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\tRA:Z:%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,rawReadBeforeBS,forQuality,Nmismatch,read_Methyl_Info,mappingStrand.c_str(), readString.c_str());
				}else 
				{
					fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,readString.c_str(),forQuality,Nmismatch,read_Methyl_Info,mappingStrand.c_str());
				}
			    }
			    funlockfile(((ARGS *)arg)->OUTFILE);
            }
		}
		if(REMOVE_DUP){
			((ARGS *)arg)->Marked_Genome[pos+G_Skip] |= Flag_rm;
			((ARGS *)arg)->Marked_GenomeE[pos+G_Skip+readString.size()] |= Flag_rm;
		}
	}
	free(s2t);
	fclose(fIS);
}

void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L) {
	char bases[] = "ACGT";

	for (int i=0; i<L; i++) {
		for (int j=0; j<4; j++) {
			for (int k=0; k<4; k++) {
				fprintf(OUTFILE_MM,"%u\t", Mismatch_Qual[i][bases[j]][bases[k]]);
			}
		} fprintf(OUTFILE_MM,"\n");
	}

}

inline unsigned char Hash(char* S)
{
	assert(false);
	unsigned char C=0;
	for (int i=2;S[i];C+=i*S[i++]);
	return C;
}

inline float calc_Entropy (string readString, int L)
{ 
	short entropy_arr[255]={0};
	//int length=strlen(readString);
	for(int i=0; i<L; i++) entropy_arr[readString[i]]++;
	
	float entropy=0.0;
	for(int i=0; i<4; i++) {
		double p = 1.0*entropy_arr["ACGT"[i]]/L;
		if(p>0) entropy-=p*log(p);
	}
  //if( entropy_arr["ACGT"[1]] + entropy_arr["ACGT"[3]] ==L ||  entropy_arr["ACGT"[0]] + entropy_arr["ACGT"[2]] == L ) entropy = 0;
	return entropy;
}

int Get_String_Length(FILE* INFILE)
{
	char Buf[BATBUF],Dummy[BATBUF],Tag[BATBUF];int L;
	fgets(Buf,BATBUF,INFILE);
	fgets(Buf,BATBUF,INFILE);
	//space-trouble in read description
	for(int i=0;Buf[i];i++) if(Buf[i]==' ') Buf[i]='_';
	sscanf(Buf,"%s%s",Dummy,Tag);	
	for (L=0;Tag[L];L++);
	rewind (INFILE);
	return L;
}
//{----------------------------------- FILE HANDLING ---------------------------------------------------------

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

//}----------------------------------- FILE HANDLING ---------------------------------------------------------


void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}
/*
string Get_RealRead(char* cig,string & read)
{
	std::string ReadR="";
	char temp[5];unsigned lens=0;
	unsigned n=0;
	while(*cig!='\0')
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("S");
			}
			lens+=atoi(temp);
			cig++;n=0;
		}else if(*cig=='M')
		{
			temp[n]='\0';std::cout<<read<<endl;
			ReadR.append(read.substr(lens,atoi(temp)));
			cig++;n=0;lens+=atoi(temp);
		}else if(*cig=='I')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("I");
			}
			lens+=atoi(temp);
			cig++;n=0;
		}else if(*cig=='D')
		{
			int i;temp[n]='\0';
			for(i=0;i<atoi(temp);i++)
			{
				ReadR.append("D");
			}
			//lens+=atoi(temp);
			cig++;n=0;
		}else
		{
			printf(" --%c-- ",*cig);
			exit(0);
		}
	}
	return ReadR;
}
*/
float Pr(float Q)
{
        assert((int)Q>=0 && (int)Q<POWLIMIT-1);
        return(1-POW10[(int)Q]);
        //printf("table: %f\tlib: %f\n",POW10[(int)Q],1-pow(10,-Q/10));
        //return (1-pow(10,-Q/10));
}
void Build_Pow10()
{
	for(int Q=0;Q<POWLIMIT;Q++)
	{
		POW10[Q]=(pow(10,-float(Q)/10));
	}
}
void fetch(const char *str, char c1, char c2, char *buf){
    const char *pl = strchr(str, c1)+1,
        *ph = strchr(str, c2);
    strncpy(buf, pl, ph-pl-1);
    buf[ph-pl] = '\0';
}
void initMem(char* temp,int size,char* temp2)
{
	temp= ( char* )malloc(size);
	strcpy(temp,temp2);
}
/*
X   for methylated C in CHG context
x   for not methylated C CHG
H   for methylated C in CHH context
h   for not methylated C in CHH context
Z   for methylated C in CpG context
z   for not methylated C in CpG context
M   for methylated C in Unknown context (CN or CHN )
U   for not methylated C in Unknown context (CN or CHN)
=    for match bases
A/T/C/G   for mismatch bases
*/
string&   replace_all(string&   str,const   string&   old_value,const   string&   new_value)   
{   
    while(true)   {   
        string::size_type   pos(0);   
        if(   (pos=str.find(old_value))!=string::npos   )   
            str.replace(pos,old_value.length(),new_value);   
        else   break;   
    }   
    return   str;   
}  

std::string m_replace(std::string str,std::string pattern,std::string dstPattern,int count)
{
    std::string retStr="";
    string::size_type pos;

    int szStr=str.length();
    int szPattern=pattern.size();
    int i=0;
    int l_count=0;
    if(-1 == count) // replace all
        count = szStr;

    for(i=0; i<szStr; i++)
    {
        pos=str.find(pattern,i);

        if(std::string::npos == pos)
            break;
        if(pos < szStr)
        {
            std::string s=str.substr(i,pos-i);
            retStr += s;
            retStr += dstPattern;
            i=pos+pattern.length()-1;
            if(++l_count >= count)
            {
                i++;
                break;
            }
        }
    }
    retStr += str.substr(i);
    return retStr;
}
int strSearch(char *str1,char *str2)
{
       int at,flag=1;
       if (strlen(str2) > strlen(str1))
       {
           at = -1;
       }
       else if (!strcmp(str1,str2))
       {
           at = 0;
       }
       else
       {
            unsigned i=0,j=0;
            for (i=0;i < strlen(str1)&&flag;)
           {
                  for (j=0;j < strlen(str2)&&flag;)
                 {
                       if (str1[i]!=str2[j])
                       {
                              i++;
                              j=0;
                       }
                       else if (str1[i]==str2[j])
                       {
                              i++;
                              j++;
                       }
                       if (j==strlen(str2))
                       {
                             flag = 0;
                       }
                       if(i==strlen(str1)) break;
                 }
            }
            at = i-j;
            if(flag==1) at=-1;
       }
       return at;
}

char* replace(char*src, char*sub, char*dst)
{
    int pos =0;
    int offset =0;
    int srcLen, subLen, dstLen;
    char*pRet = NULL;

    srcLen = strlen(src);
    subLen = strlen(sub);
    dstLen = strlen(dst);
    pRet = (char*)malloc(srcLen + dstLen - subLen +1);//(memory)if (NULL != pRet)
    {
        pos = strstr(src, sub) - src;
        memcpy(pRet, src, pos);
        offset += pos;
        memcpy(pRet + offset, dst, dstLen);
        offset += dstLen;
        memcpy(pRet + offset, src + pos + subLen, srcLen - pos - subLen);
        offset += srcLen - pos - subLen;
        *(pRet + offset) ='\0';
    }
    return pRet;
}
int Get_ED(std::string & S)//edit distance
{
	int ed=0;
	int i;
	for(i=0;i<S.size();i++)//momomo
	{
		if(S[i]=='I' || S[i]=='D') ed++;
	}
	return ed;
}
// 
void MinandSec(unsigned a[MAX_HITS_ALLOWED][10], int left, int right, int & min, int & second,int & loci)
{
	if(left == right)
	{
		min = a[left][4] ;
		second =  INT_MIN;
		loci=left;
	}
	else if(left +1== right)
	{
		min = a[left][4] < a[right][4] ? a[left][4] : a[right][4] ;
		second = a[left][4] > a[right][4] ? a[left][4] : a[right][4] ;
		loci=a[left][4] < a[right][4] ? left : right ;
	}
	else
	{
		int mid =(right + left) /2 ;

		int leftmin ;
		int leftsecond ;
		int leftloci;
		MinandSec(a, left, mid, leftmin, leftsecond,leftloci) ;

		int rightmin ;
		int rightsecond ;
		int rightloci;
		MinandSec(a, mid +1, right, rightmin, rightsecond,rightloci) ;

		if (leftmin < rightmin)
		{
			min = leftmin ;
			second = leftsecond < rightmin ? leftsecond : rightmin ;
			loci=leftloci;
		}
		else
		{
			min = rightmin ;
			second = leftmin > rightsecond ? rightsecond : leftmin ;
			loci=rightloci;
		}
	}
	return;
}
/*
void MaxandSec(int a[], int left, int right, int&max, int&second)
{
	if(left == right)
	{
		max = a[left] ;
		second =  INT_MIN;
	}
	else if(left +1== right)
	{
		max = a[left] > a[right] ? a[left] : a[right] ;
		second = a[left] < a[right] ? a[left] : a[right] ;
	}
	else
	{
		int mid =(right + left) /2 ;

		int leftmax ;
		int leftsecond ;
		MaxandSec(a, left, mid, leftmax, leftsecond) ;

		int rightmax ;
		int rightsecond ;
		MaxandSec(a, mid +1, right, rightmax, rightsecond) ;

		if (leftmax > rightmax)
		{
			max = leftmax ;
			second = leftsecond > rightmax ? leftsecond : rightmax ;
		}
		else
		{
			max = rightmax ;
			second = leftmax < rightsecond ? rightsecond : leftmax ;
		}
	}
}
*/

string remove_soft_split(string & cigar,int & Read_Len,int & pos)
{

	const char* cig=cigar.c_str();
	char temp[20];
	int n=0,lens=0,length=0,RLens=0,Glens=0;
	int genome_move_size=0; int headClip=0;int moveSize=0;
	char cigar_rm[1000];*cigar_rm=0;
	char buffer_cigar[100];*buffer_cigar=0; bool cigar_remove=false;
	char upper='N';char upper_cigar[1024]="\0";int upper_length=0;
	while(*cig!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			if(length==0) cigar_remove=true;
			if(length>0) 
			{
				if(upper=='M')
					sprintf(upper_cigar,"%dM",length+upper_length);
				else
					sprintf(upper_cigar,"%dM",length);
				
				if(upper=='N')
					pos-=length;
			}

			lens+=length;
            RLens+=length;
            cig++;n=0;
            upper='M';upper_length=length;
		}else if(*cig=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			if(length==0) cigar_remove=true;
			if(length>0) 
			{
				if(upper=='M')
					sprintf(upper_cigar,"%dM",length+upper_length);
				else
					sprintf(upper_cigar,"%dM",length);
			}

			cig++;n=0;lens+=length;Glens+=length;RLens+=length;
			char buf[1024]="\0";
			sprintf( buf , "%dM",Read_Len);
			char buf_tmp[1024]="\0";
			sprintf( buf_tmp, "%dM",length);
			if( !strcmp(buf,  buf_tmp ) ) 
			{ 
				if(!strcmp(buf, cig))
					break;
				else
				{
					cigar_remove=true;
					strcpy(cigar_rm, buffer_cigar);
				}
				
			}
			upper='M';upper_length=length;
		}else if(*cig=='I')
		{
			if(upper=='M')
				strcat(cigar_rm, upper_cigar);
			int i;temp[n]='\0';int length=atoi(temp);

			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1I");
				else if(length==2) sprintf(buffer_cigar,"2I");
				else sprintf(buffer_cigar,"%dI",length);
				strcat(cigar_rm, buffer_cigar);
			}

			lens+=length;RLens+=length;
			cig++;n=0;genome_move_size-=length;moveSize++;
			upper='I';
		}else if(*cig=='D')
		{
			if(upper=='M')
				strcat(cigar_rm, upper_cigar);
			int i;temp[n]='\0';int length=atoi(temp);

			if(length>0) 
			{
				if(length==1) sprintf(buffer_cigar,"1D");
				else if(length==2) sprintf(buffer_cigar,"2D");
				else sprintf(buffer_cigar,"%dD",length);
				strcat(cigar_rm, buffer_cigar);
			}

			Glens+=length;RLens+=length;
			cig++;n=0;genome_move_size+=length;moveSize++;
			upper='D';
		}else
		{
			printf(" --%d%c-- ",atoi(temp),*cig);
			break;
		}
	}
	if(upper=='M')
		strcat(cigar_rm, upper_cigar);
	string s(cigar_rm);
	return s;
}

void ReverseC_Context(char* Dest,const char* seq,int & stringlength)
{
	if(stringlength!=5 || strlen(seq)!=5) {
		fprintf(stderr, "\nError string %d %d %s\n", stringlength, strlen(seq), seq);
		exit(0);
	}
	
        for (int i=stringlength-1;i>=0;i--)
        {
                *Dest=Char2Comp[seq[i]];
                Dest++;
        }
        *Dest=0;
        //return Dest;
}

int fprintf_time(FILE *stream, const char *format, ...)
{
        time_t timer;
        char buffer[26];
        struct tm* tm_info;
        va_list arg;
        int done;

        time(&timer);
        tm_info = localtime(&timer);
        strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", tm_info);

        fprintf(stream, "[DM::calmeth] [%s] ", buffer);

        va_start(arg, format);
        done = vfprintf(stream, format, arg);
        va_end(arg);

        return done;

}
