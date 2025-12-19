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
#include <unordered_map>
#include <vector>
#include <sys/stat.h>
#include <stdint.h>
#include <inttypes.h>
#include <cerrno>
#include <thread>
#include <mutex>
#include <atomic>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <dirent.h>
extern "C" {
#include "dmCommon.h"
}

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

//#define DEBUG 2
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
   FILE* OUTFILEMS;
   FILE* methOutFp;
   FILE* samINFILE;
   FILE* bedFILE;
   samFile* BamInFile;
   bam1_t *b;
   bam_hdr_t *header;
   int ThreadID;
   off64_t File_Size;
   char* processChr;
   char* INbamfilename;
   char methOutPath[PATH_MAX];
   uint64_t methOutLines;
} ARGS;
static void freeMethArrays(ARGS &args);
bool RELESEM = true;// false
bool printheader = true;
using namespace std;
bool Collision=false;
map <string,int> String_Hash;
float ENTROPY_CUTOFF=0;
int skipOverlap = 1; //
int countreadC = 0;
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
void Reverse_Comp(char* Dest, char* seq);
void onlyComp(char* Dest, char* seq);
//---------------------------------------------------------------------------------------------------------------

static size_t estimateAlignedColumns(const std::string &cigar) {
    size_t readBases = 0;
    size_t insertBases = 0;
    size_t deleteBases = 0;
    size_t colBases = 0;
    char numbuf[32];
    int n = 0;
    for(const char *c = cigar.c_str(); *c; ++c) {
        if(isdigit(*c)) {
            if(n < (int)sizeof(numbuf) - 1) {
                numbuf[n++] = *c;
            } else {
                return cigar.size() + 1024; // fallback for pathological inputs
            }
            continue;
        }
        numbuf[n] = '\0';
        n = 0;
        int len = atoi(numbuf);
        switch(*c) {
            case 'M': case '=': case 'X':
                readBases += len;
                colBases += len;
                break;
            case 'I':
                readBases += len;
                insertBases += len;
                colBases += len;
                break;
            case 'D': case 'N':
                deleteBases += len;
                colBases += len;
                break;
            case 'S': case 'H':
                readBases += len;
                break;
            default:
                break;
        }
    }
    size_t aligned = readBases + deleteBases;
    if(aligned < colBases) aligned = colBases;
    aligned += insertBases;
    return aligned + 16; // safety padding
}
struct BinTask {
    std::string chrom;
    int64_t start;
    int64_t end;
    size_t index;
};

struct BinTaskStats {
    BinTask task;
    uint64_t reads = 0;
    uint64_t records = 0;
    int status = 0;
};

static bool gDebugMode = false;
static uint64_t gDmRecordsWritten = 0;
static uint64_t gDmRecordsWrittenGch = 0;
static uint64_t gGenomeSizeBases = 0;
static bool gDmWritersClosed = false;
static std::mutex gDmWriterMutex;
static std::once_flag gMethArraysCleanupFlag;

struct FilterStats {
    std::atomic<uint64_t> totalReads{0};
    std::atomic<uint64_t> processReadAccepted{0};
    std::atomic<uint64_t> bamFiltered{0};
    std::atomic<uint64_t> mapqFiltered{0};
    std::atomic<uint64_t> hashMiss{0};
    std::atomic<uint64_t> mismatchFiltered{0};
    std::atomic<uint64_t> refOob{0};
    std::atomic<uint64_t> refNonACGT{0};
    std::atomic<uint64_t> enterMethyl{0};
    std::atomic<uint64_t> methylCalls{0};
    std::atomic<uint64_t> iterNull{0};
};

static FilterStats gFilterStats;

static void maybeLogFilterStats() {
    static std::atomic<uint64_t> lastLogged{0};
    const uint64_t cur = gFilterStats.totalReads.load();
    if(cur - lastLogged.load() >= 1000000) {
        lastLogged.store(cur);
        fprintf(stderr,
                "[filter] reads=%" PRIu64 " accepted=%" PRIu64 " bamFiltered=%" PRIu64 " mapq=%" PRIu64
                " hashMiss=%" PRIu64 " mismatch=%" PRIu64 " refOOB=%" PRIu64 " refNonACGT=%" PRIu64
                " enterMethyl=%" PRIu64 " methylCalls=%" PRIu64 " iterNull=%" PRIu64 "\n",
                cur,
                gFilterStats.processReadAccepted.load(),
                gFilterStats.bamFiltered.load(),
                gFilterStats.mapqFiltered.load(),
                gFilterStats.hashMiss.load(),
                gFilterStats.mismatchFiltered.load(),
                gFilterStats.refOob.load(),
                gFilterStats.refNonACGT.load(),
                gFilterStats.enterMethyl.load(),
                gFilterStats.methylCalls.load(),
                gFilterStats.iterNull.load());
    }
}

struct BinPartFileHeader {
    uint32_t tid;
    uint64_t start;
    uint64_t end;
    uint32_t nRecords;
};

struct BinPartLocator {
    int shard;
    uint64_t offset;
    uint32_t nRecords;
    uint32_t tid;
    uint64_t start;
    uint64_t end;
};

struct BinPartRecord {
    uint32_t pos;
    uint32_t end;
    float value;
    uint16_t coverage;
    uint8_t strand;
    uint8_t context;
};

static int runBinChunkDebug(const std::string &bamPath, const std::vector<BinTask> &tasks, int threads, bool debugMode);
static int materializeBinTasksToBed(const std::string &genomePath, const std::string &targetChrom, int binSize,
                                    std::string &bedOut, size_t &taskCount, bool debugMode);
static int loadBinTasksFromBed(const std::string &bedPath, std::vector<BinTask> &tasks, bool debugMode);
static int validateDmFile(const std::string &dmPath, bool verbose);
static int runBinTasksToParts(const std::string &bamPath, const std::vector<BinTask> &tasks, int threads,
                              const std::string &partDir, bool debugMode, std::vector<BinPartLocator> &locators);
static int mergeBinPartsToDm(const std::string &genomePath, const std::vector<BinTask> &tasks,
                             const std::vector<BinPartLocator> &locators, int shardCount,
                             const std::string &partDir, const std::string &dmPath, int zoomlevel,
                             bool debugMode, bool validateOutput);

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
int QualCut=30;
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
int printtxt=0;
int Mcoverage=4;
int maxcoverage=1000;
int binspan=50000;
int nCs=1;
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
string hsPrefix="hseq.txt";
#define MAX_LINE_PRINT 1000000
#define BUFSIZE 1000000
#define MAX_CHROM 1000000
binaMethFile_t *fp;
binaMethFile_t *fp_gch; //GCH sites (GCA/GCT/GCC) were used to analyze chromatin accessibility
string contextfilter="C";
string tech="BS-Seq";
std::unordered_map<std::string, int> myMap;
int onlyPHead = 0;
int PHead = 0;
int rrbs=0;
int printmethstate = 0;
int onlyM = 0;

static void destroyOverlapBlock(bmOverlapBlock_t *block) {
    if(!block) return;
    if(block->size) free(block->size);
    if(block->offset) free(block->offset);
    free(block);
}

static void closeDmOutputs() {
    std::lock_guard<std::mutex> lock(gDmWriterMutex);
    if(gDmWritersClosed) return;
    gDmWritersClosed = true;
    if(gDebugMode) {
        fprintf(stderr, "[dm-writer] closing outputs fp=%p fp_gch=%p\n", (void*)fp, (void*)fp_gch);
    }
    if(fp && fp->writeBuffer) {
        if(gDebugMode) fprintf(stderr, "[dm-writer] closing main dm writer\n");
        bmClose(fp);
        fp = NULL;
    }
    if(fp_gch && fp_gch->writeBuffer) {
        if(gDebugMode) fprintf(stderr, "[dm-writer] closing GCH dm writer\n");
        bmClose(fp_gch);
        fp_gch = NULL;
    }
}

static int runBinChunkDebug(const std::string &bamPath, const std::vector<BinTask> &tasks, int threads, bool debugMode) {
    if(bamPath.empty()) {
        fprintf(stderr, "[bin] input BAM/SAM is required for --chunk-by bin mode\n");
        return 1;
    }
    if(tasks.empty()) {
        fprintf(stderr, "[bin] no bin tasks to process\n");
        return 1;
    }

    if(threads < 1) threads = 1;
    if(debugMode) {
        fprintf(stderr, "[bin] initializing %zu tasks across %d thread(s)\n", tasks.size(), threads);
    }

    std::atomic<size_t> nextTask(0);
    std::vector<BinTaskStats> stats(tasks.size());
    std::mutex logMutex;
    std::atomic<int> workerError(0);

    auto worker = [&](int workerId) {
        samFile *bam = sam_open(bamPath.c_str(), "rb");
        if(!bam) {
            fprintf(stderr, "[bin-worker-%d] failed to open %s\n", workerId, bamPath.c_str());
            workerError.store(1);
            return;
        }
        bam_hdr_t *header = sam_hdr_read(bam);
        if(!header) {
            fprintf(stderr, "[bin-worker-%d] failed to read BAM header from %s\n", workerId, bamPath.c_str());
            sam_close(bam);
            workerError.store(1);
            return;
        }
        hts_idx_t *idx = sam_index_load(bam, bamPath.c_str());
        if(!idx) {
            fprintf(stderr, "[bin-worker-%d] failed to load BAM index for %s\n", workerId, bamPath.c_str());
            bam_hdr_destroy(header);
            sam_close(bam);
            workerError.store(1);
            return;
        }

        bam1_t *b = bam_init1();
        if(!b) {
            fprintf(stderr, "[bin-worker-%d] failed to allocate bam1_t\n", workerId);
            hts_idx_destroy(idx);
            bam_hdr_destroy(header);
            sam_close(bam);
            workerError.store(1);
            return;
        }

        while(true) {
            size_t taskId = nextTask.fetch_add(1);
            if(taskId >= tasks.size()) break;
            const BinTask &task = tasks[taskId];
            BinTaskStats result;
            result.task = task;

            int tid = bam_name2id(header, task.chrom.c_str());
            if(tid < 0) {
                result.status = 1;
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] chrom %s not found in BAM header\n", workerId, task.chrom.c_str());
                stats[taskId] = result;
                continue;
            }

            hts_itr_t *iter = sam_itr_queryi(idx, tid, task.start, task.end);
            if(!iter) {
                result.status = 1;
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] failed to create iterator for %s:%" PRId64 "-%" PRId64 "\n", workerId,
                        task.chrom.c_str(), task.start, task.end);
                stats[taskId] = result;
                continue;
            }

            {
                std::lock_guard<std::mutex> lock(logMutex);
                if(debugMode) {
                    fprintf(stderr, "[bin-worker-%d] start task %zu %s:%" PRId64 "-%" PRId64 "\n", workerId, task.index,
                            task.chrom.c_str(), task.start, task.end);
                }
            }

            while(sam_itr_next(bam, iter, b) >= 0) {
                result.reads++;
                result.records++;
            }

            hts_itr_destroy(iter);

            {
                std::lock_guard<std::mutex> lock(logMutex);
                if(debugMode) {
                    fprintf(stderr, "[bin-worker-%d] done task %zu %s:%" PRId64 "-%" PRId64 " reads=%" PRIu64 "\n", workerId,
                            task.index, task.chrom.c_str(), task.start, task.end, result.reads);
                }
            }

            stats[taskId] = result;
        }

        bam_destroy1(b);
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam);
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for(int t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for(auto &th : pool) {
        th.join();
    }

    if(workerError.load() != 0) {
        fprintf(stderr, "[bin] worker error detected; aborting\n");
        return 1;
    }

    for(size_t i = 0; i < stats.size(); ++i) {
        if(stats[i].status != 0) {
            fprintf(stderr, "[bin] task %zu (%s:%" PRId64 "-%" PRId64 ") failed\n", i, stats[i].task.chrom.c_str(),
                    stats[i].task.start, stats[i].task.end);
            return 1;
        }
    }

    if(debugMode) {
        uint64_t totalReads = 0;
        for(const auto &s : stats) totalReads += s.reads;
        fprintf(stderr, "[bin] completed %zu tasks; total reads scanned=%" PRIu64 "\n", stats.size(), totalReads);
    }

    return 0;
}

static int materializeBinTasksToBed(const std::string &genomePath, const std::string &targetChrom, int binSize,
                                    std::string &bedOut, size_t &taskCount, bool debugMode) {
    taskCount = 0;
    if(genomePath.empty()) {
        fprintf(stderr, "[bin] --genome is required for --chunk-by bin mode\n");
        return 1;
    }
    if(binSize <= 0) {
        fprintf(stderr, "[bin] bin-size must be positive\n");
        return 1;
    }

    char tmpl[] = "/tmp/bam2dm_binsXXXXXX";
    int fd = mkstemp(tmpl);
    if(fd == -1) {
        fprintf(stderr, "[bin] failed to create temporary bed for bin tasks: %s\n", strerror(errno));
        return 1;
    }
    FILE *out = fdopen(fd, "w+");
    if(!out) {
        fprintf(stderr, "[bin] failed to open temp bed stream: %s\n", strerror(errno));
        close(fd);
        return 1;
    }

    faidx_t *fai = fai_load(genomePath.c_str());
    if(!fai) {
        fprintf(stderr, "[bin] failed to load genome index: %s\n", genomePath.c_str());
        fclose(out);
        unlink(tmpl);
        return 1;
    }

    const int nseq = faidx_nseq(fai);
    for(int i = 0; i < nseq; ++i) {
        const char *name = faidx_iseq(fai, i);
        if(!name) continue;
        if(strcmp(targetChrom.c_str(), "NAN-mm") != 0 && targetChrom != name) {
            continue;
        }
        const int64_t len = faidx_seq_len(fai, name);
        if(len <= 0) continue;
        for(int64_t start = 0; start < len; start += binSize) {
            int64_t end = start + binSize;
            if(end > len) end = len;
            fprintf(out, "%s\t%" PRId64 "\t%" PRId64 "\n", name, start, end);
            ++taskCount;
        }
    }

    fai_destroy(fai);
    fflush(out);
    fseek(out, 0, SEEK_SET);
    bedOut.assign(tmpl);

    if(debugMode) {
        fprintf(stderr, "[bin] materialized %zu bin task(s) to %s\n", taskCount, bedOut.c_str());
    }

    fclose(out);
    return 0;
}

static int loadBinTasksFromBed(const std::string &bedPath, std::vector<BinTask> &tasks, bool debugMode) {
    tasks.clear();
    std::ifstream in(bedPath.c_str());
    if(!in.is_open()) {
        fprintf(stderr, "[bin] failed to open bed task file %s\n", bedPath.c_str());
        return 1;
    }
    std::string line;
    size_t idx = 0;
    while(std::getline(in, line)) {
        if(line.empty()) continue;
        std::istringstream iss(line);
        std::string chrom;
        int64_t start = 0, end = 0;
        if(!(iss >> chrom >> start >> end)) {
            fprintf(stderr, "[bin] malformed bed line: %s\n", line.c_str());
            return 1;
        }
        if(start < 0 || end < 0 || end < start) {
            fprintf(stderr, "[bin] invalid interval in bed: %s\n", line.c_str());
            return 1;
        }
        BinTask task;
        task.chrom = chrom;
        task.start = start;
        task.end = end;
        task.index = idx++;
        tasks.push_back(task);
    }
    if(debugMode) {
        fprintf(stderr, "[bin] loaded %zu tasks from %s\n", tasks.size(), bedPath.c_str());
    }
    return tasks.empty() ? 1 : 0;
}

static int runBinTasksToParts(const std::string &bamPath, const std::vector<BinTask> &tasks, int threads,
                              const std::string &partDir, bool debugMode, std::vector<BinPartLocator> &locators) {
    if(bamPath.empty()) {
        fprintf(stderr, "[bin] input BAM/SAM is required for --chunk-by bin mode\n");
        return 1;
    }
    if(tasks.empty()) {
        fprintf(stderr, "[bin] no bin tasks to process\n");
        return 1;
    }

    if(threads < 1) threads = 1;
    if(debugMode) {
        fprintf(stderr, "[bin] initializing %zu tasks across %d thread(s) (shards=%d dir=%s)\n", tasks.size(), threads, threads,
                partDir.c_str());
    }

    if(partDir.empty()) {
        fprintf(stderr, "[bin] part directory is empty\n");
        return 1;
    }
    if(mkdir(partDir.c_str(), 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "[bin] failed to create part directory %s: %s\n", partDir.c_str(), strerror(errno));
        return 1;
    }

    std::atomic<size_t> nextTask(0);
    std::atomic<int> workerError(0);
    std::mutex logMutex;
    locators.assign(tasks.size(), {-1, 0, 0, 0, 0, 0});

    std::vector<std::string> shardPaths;
    shardPaths.reserve(threads);
    for(int t = 0; t < threads; ++t) {
        char tmpPath[PATH_MAX];
        snprintf(tmpPath, sizeof(tmpPath), "%s/thread_%d.tmp", partDir.c_str(), t);
        shardPaths.emplace_back(tmpPath);
    }

    auto worker = [&](int workerId) {
        FILE *shard = fopen(shardPaths[workerId].c_str(), "wb+");
        if(!shard) {
            std::lock_guard<std::mutex> lock(logMutex);
            fprintf(stderr, "[bin-worker-%d] failed to open shard %s: %s\n", workerId, shardPaths[workerId].c_str(),
                    strerror(errno));
            workerError.store(1);
            return;
        }
        samFile *bam = sam_open(bamPath.c_str(), "rb");
        if(!bam) {
            fprintf(stderr, "[bin-worker-%d] failed to open %s\n", workerId, bamPath.c_str());
            workerError.store(1);
            fclose(shard);
            return;
        }
        bam_hdr_t *header = sam_hdr_read(bam);
        if(!header) {
            fprintf(stderr, "[bin-worker-%d] failed to read BAM header from %s\n", workerId, bamPath.c_str());
            sam_close(bam);
            workerError.store(1);
            fclose(shard);
            return;
        }
        hts_idx_t *idx = sam_index_load(bam, bamPath.c_str());
        if(!idx) {
            fprintf(stderr, "[bin-worker-%d] failed to load BAM index for %s\n", workerId, bamPath.c_str());
            bam_hdr_destroy(header);
            sam_close(bam);
            workerError.store(1);
            fclose(shard);
            return;
        }

        bam1_t *b = bam_init1();
        if(!b) {
            fprintf(stderr, "[bin-worker-%d] failed to allocate bam1_t\n", workerId);
            hts_idx_destroy(idx);
            bam_hdr_destroy(header);
            sam_close(bam);
            fclose(shard);
            workerError.store(1);
            return;
        }

        while(true) {
            size_t taskId = nextTask.fetch_add(1);
            if(taskId >= tasks.size()) break;
            const BinTask &task = tasks[taskId];
            int tid = bam_name2id(header, task.chrom.c_str());
            if(tid < 0) {
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] chrom %s not found in BAM header\n", workerId, task.chrom.c_str());
                workerError.store(1);
                continue;
            }

            hts_itr_t *iter = sam_itr_queryi(idx, tid, task.start, task.end);
            if(!iter) {
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] failed to create iterator for %s:%" PRId64 "-%" PRId64 "\n", workerId,
                        task.chrom.c_str(), task.start, task.end);
                workerError.store(1);
                continue;
            }

            std::vector<BinPartRecord> records;
            while(sam_itr_next(bam, iter, b) >= 0) {
                const int64_t pos = b->core.pos;
                if(pos < task.start || pos >= task.end) continue;
                BinPartRecord rec{};
                rec.pos = static_cast<uint32_t>(pos);
                rec.end = rec.pos + 1;
                rec.value = 1.0f;
                rec.coverage = 1;
                rec.strand = (b->core.flag & BAM_FREVERSE) ? 1 : 0;
                rec.context = 0;
                records.push_back(rec);
            }
            hts_itr_destroy(iter);

            std::sort(records.begin(), records.end(), [](const BinPartRecord &a, const BinPartRecord &b) {
                if(a.pos == b.pos) return a.end < b.end;
                return a.pos < b.pos;
            });

            BinPartLocator loc{};
            loc.shard = workerId;
            loc.tid = static_cast<uint32_t>(tid);
            loc.start = static_cast<uint64_t>(task.start);
            loc.end = static_cast<uint64_t>(task.end);
            loc.nRecords = static_cast<uint32_t>(records.size());

            off64_t off = ftello64(shard);
            if(off < 0) {
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] failed to get offset for shard %s\n", workerId, shardPaths[workerId].c_str());
                workerError.store(1);
            } else {
                loc.offset = static_cast<uint64_t>(off);
                BinPartFileHeader hdr{};
                hdr.tid = loc.tid;
                hdr.start = loc.start;
                hdr.end = loc.end;
                hdr.nRecords = loc.nRecords;
                size_t wrote = fwrite(&hdr, sizeof(BinPartFileHeader), 1, shard);
                if(wrote != 1) {
                    std::lock_guard<std::mutex> lock(logMutex);
                    fprintf(stderr, "[bin-worker-%d] failed to write header for shard %s\n", workerId, shardPaths[workerId].c_str());
                    workerError.store(1);
                } else if(!records.empty()) {
                    wrote = fwrite(records.data(), sizeof(BinPartRecord), records.size(), shard);
                    if(wrote != records.size()) {
                        std::lock_guard<std::mutex> lock(logMutex);
                        fprintf(stderr, "[bin-worker-%d] failed to write records for shard %s\n", workerId, shardPaths[workerId].c_str());
                        workerError.store(1);
                    }
                }
                locators[task.index] = loc;
            }

            if(debugMode) {
                std::lock_guard<std::mutex> lock(logMutex);
                fprintf(stderr, "[bin-worker-%d] task %zu %s:%" PRId64 "-%" PRId64 " records=%zu\n", workerId, task.index,
                        task.chrom.c_str(), task.start, task.end, records.size());
            }
        }

        bam_destroy1(b);
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam);
        fclose(shard);
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for(int t = 0; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for(auto &th : pool) th.join();

    return workerError.load();
}

static int mergeBinPartsToDm(const std::string &genomePath, const std::vector<BinTask> &tasks,
                             const std::vector<BinPartLocator> &locators, int shardCount,
                             const std::string &partDir, const std::string &dmPath, int zoomlevel,
                             bool debugMode, bool validateOutput) {
    if(tasks.empty()) {
        fprintf(stderr, "[bin] no tasks to merge\n");
        return 1;
    }

    if(locators.size() != tasks.size()) {
        fprintf(stderr, "[bin] locator map does not match task list\n");
        return 1;
    }

    faidx_t *fai = fai_load(genomePath.c_str());
    if(!fai) {
        fprintf(stderr, "[bin] failed to load genome index %s\n", genomePath.c_str());
        return 1;
    }
    int nseq = faidx_nseq(fai);
    if(nseq <= 0) {
        fprintf(stderr, "[bin] empty genome index %s\n", genomePath.c_str());
        fai_destroy(fai);
        return 1;
    }

    char **chroms = (char **)calloc(nseq, sizeof(char *));
    uint32_t *lens = (uint32_t *)calloc(nseq, sizeof(uint32_t));
    if(!chroms || !lens) {
        fprintf(stderr, "[bin] failed to allocate chrom list for dm writer\n");
        fai_destroy(fai);
        free(chroms); free(lens);
        return 1;
    }
    for(int i = 0; i < nseq; ++i) {
        const char *name = faidx_iseq(fai, i);
        chroms[i] = strdup(name);
        lens[i] = static_cast<uint32_t>(faidx_seq_len(fai, name));
    }

    binaMethFile_t *out = (binaMethFile_t*)bmOpen((char*)dmPath.c_str(), NULL, "w");
    if(!out) {
        fprintf(stderr, "[bin] failed to open output dm %s\n", dmPath.c_str());
        for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
        free(chroms); free(lens);
        fai_destroy(fai);
        return 1;
    }
    if(bmCreateHdr(out, zoomlevel)) {
        fprintf(stderr, "[bin] failed to create dm header for %s\n", dmPath.c_str());
        bmClose(out);
        for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
        free(chroms); free(lens);
        fai_destroy(fai);
        return 1;
    }
    out->cl = bmCreateChromList(chroms, lens, nseq);
    if(!out->cl) {
        fprintf(stderr, "[bin] failed to create chrom list for %s\n", dmPath.c_str());
        bmClose(out);
        for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
        free(chroms); free(lens);
        fai_destroy(fai);
        return 1;
    }
    if(bmWriteHdr(out)) {
        fprintf(stderr, "[bin] failed to write header for %s\n", dmPath.c_str());
        bmClose(out);
        for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
        free(chroms); free(lens);
        fai_destroy(fai);
        return 1;
    }

    std::vector<BinTask> ordered = tasks;
    std::sort(ordered.begin(), ordered.end(), [](const BinTask &a, const BinTask &b) {
        if(a.chrom == b.chrom) return a.start < b.start;
        return a.chrom < b.chrom;
    });

    std::vector<FILE*> shardReaders(shardCount, NULL);
    for(int s = 0; s < shardCount; ++s) {
        char shardPath[PATH_MAX];
        snprintf(shardPath, sizeof(shardPath), "%s/thread_%d.tmp", partDir.c_str(), s);
        shardReaders[s] = fopen(shardPath, "rb");
        if(!shardReaders[s]) {
            fprintf(stderr, "[bin] failed to open shard %s\n", shardPath);
            for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
            free(chroms); free(lens);
            fai_destroy(fai);
            bmClose(out);
            return 1;
        }
    }

    std::vector<char*> chromBuf(MAX_LINE_PRINT, NULL);
    std::vector<uint32_t> startBuf(MAX_LINE_PRINT, 0);
    std::vector<uint32_t> endBuf(MAX_LINE_PRINT, 0);
    std::vector<float> valueBuf(MAX_LINE_PRINT, 0.0f);
    std::vector<uint16_t> covBuf(MAX_LINE_PRINT, 0);
    std::vector<uint8_t> strandBuf(MAX_LINE_PRINT, 0);
    std::vector<uint8_t> contextBuf(MAX_LINE_PRINT, 0);
    std::vector<char*> entryBuf(MAX_LINE_PRINT, NULL);

    for(const auto &task : ordered) {
        const BinPartLocator &loc = locators[task.index];
        if(loc.shard < 0 || loc.shard >= shardCount) {
            fprintf(stderr, "[bin] missing shard assignment for task %zu (%s:%" PRId64 "-%" PRId64 ")\n",
                    task.index, task.chrom.c_str(), task.start, task.end);
            bmClose(out);
            for(auto fh : shardReaders) if(fh) fclose(fh);
            for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
            free(chroms); free(lens);
            fai_destroy(fai);
            return 1;
        }
        FILE *in = shardReaders[loc.shard];
        if(fseeko64(in, static_cast<off64_t>(loc.offset), SEEK_SET) != 0) {
            fprintf(stderr, "[bin] failed to seek shard %d for task %zu\n", loc.shard, task.index);
            bmClose(out);
            for(auto fh : shardReaders) if(fh) fclose(fh);
            for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
            free(chroms); free(lens);
            fai_destroy(fai);
            return 1;
        }
        BinPartFileHeader hdr{};
        if(fread(&hdr, sizeof(BinPartFileHeader), 1, in) != 1) {
            fprintf(stderr, "[bin] failed to read header for task %zu from shard %d\n", task.index, loc.shard);
            bmClose(out);
            for(auto fh : shardReaders) if(fh) fclose(fh);
            for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
            free(chroms); free(lens);
            fai_destroy(fai);
            return 1;
        }
        if(hdr.nRecords != loc.nRecords || hdr.start != loc.start || hdr.end != loc.end || hdr.tid != loc.tid) {
            fprintf(stderr, "[bin] shard metadata mismatch for task %zu (shard %d)\n", task.index, loc.shard);
            bmClose(out);
            for(auto fh : shardReaders) if(fh) fclose(fh);
            for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
            free(chroms); free(lens);
            fai_destroy(fai);
            return 1;
        }

        if(hdr.nRecords > 0) {
            std::vector<BinPartRecord> recs(hdr.nRecords);
            size_t rd = fread(recs.data(), sizeof(BinPartRecord), hdr.nRecords, in);
            if(rd != hdr.nRecords) {
                fprintf(stderr, "[bin] truncated part file for task %zu on shard %d\n", task.index, loc.shard);
                bmClose(out);
                for(auto fh : shardReaders) if(fh) fclose(fh);
                for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
                free(chroms); free(lens);
                fai_destroy(fai);
                return 1;
            }

            size_t offset = 0;
            while(offset < recs.size()) {
                size_t chunk = std::min(static_cast<size_t>(MAX_LINE_PRINT), recs.size() - offset);
                for(size_t i = 0; i < chunk; ++i) {
                    const BinPartRecord &r = recs[offset + i];
                    chromBuf[i] = chroms[hdr.tid];
                    startBuf[i] = r.pos;
                    endBuf[i] = r.end;
                    valueBuf[i] = r.value;
                    covBuf[i] = r.coverage;
                    strandBuf[i] = r.strand;
                    contextBuf[i] = r.context;
                }
                int response = bmAddIntervals(out, chromBuf.data(), startBuf.data(), endBuf.data(), valueBuf.data(),
                                              covBuf.data(), strandBuf.data(), contextBuf.data(), entryBuf.data(),
                                              static_cast<uint32_t>(chunk));
                if(response != 0) {
                    fprintf(stderr, "[bin] bmAddIntervals failed for %s (code %d)\n", chroms[hdr.tid], response);
                    fclose(in);
                    bmClose(out);
                    for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
                    free(chroms); free(lens);
                    fai_destroy(fai);
                    return 1;
                }
                offset += chunk;
            }
        }
    }

    for(auto fh : shardReaders) if(fh) fclose(fh);
    bmClose(out);
    for(int i = 0; i < nseq; ++i) if(chroms[i]) free(chroms[i]);
    free(chroms); free(lens);
    fai_destroy(fai);

    if(validateOutput) {
        return validateDmFile(dmPath, debugMode);
    }
    return 0;
}

static int validateDmFile(const std::string &dmPath, bool verbose) {
    struct stat st{};
    if(stat(dmPath.c_str(), &st) != 0) {
        fprintf(stderr, "[dmcheck] failed to stat %s: %s\n", dmPath.c_str(), strerror(errno));
        return 1;
    }
    const uint64_t fileSize = static_cast<uint64_t>(st.st_size);
    binaMethFile_t *fp = bmOpen((char*)dmPath.c_str(), NULL, "r");
    if(!fp || !fp->hdr) {
        fprintf(stderr, "[dmcheck] unable to open %s as dm\n", dmPath.c_str());
        if(fp) bmClose(fp);
        return 2;
    }
    if(!fp->cl || fp->cl->nKeys == 0) {
        fprintf(stderr, "[dmcheck] missing chromosome list in %s (nKeys=%u)\n", dmPath.c_str(), fp->cl ? fp->cl->nKeys : 0);
        bmClose(fp);
        return 2;
    }

    uint32_t magic = 0;
    if(bmSetPos(fp, 0) || bmRead(&magic, sizeof(uint32_t), 1, fp) != 1) {
        fprintf(stderr, "[dmcheck] failed to read magic from %s\n", dmPath.c_str());
        bmClose(fp);
        return 3;
    }
    if(magic != BIGWIG_MAGIC) {
        fprintf(stderr, "[dmcheck] magic mismatch in %s (got 0x%x)\n", dmPath.c_str(), magic);
        bmClose(fp);
        return 4;
    }
    if(!(fp->hdr->version & BM_MAGIC)) {
        fprintf(stderr, "[dmcheck] version/type missing BM_MAGIC bit in %s\n", dmPath.c_str());
        bmClose(fp);
        return 5;
    }
    if(fp->hdr->dataOffset >= fileSize || fp->hdr->indexOffset >= fileSize) {
        fprintf(stderr, "[dmcheck] header offsets exceed file size (%s)\n", dmPath.c_str());
        bmClose(fp);
        return 6;
    }
    if(fp->hdr->ctOffset >= fileSize || (fp->hdr->summaryOffset && fp->hdr->summaryOffset + sizeof(uint64_t) * 2 > fileSize)) {
        fprintf(stderr, "[dmcheck] chromosome or summary offset outside file for %s\n", dmPath.c_str());
        bmClose(fp);
        return 7;
    }
    if(!fp->idx || !fp->idx->root) {
        fprintf(stderr, "[dmcheck] index missing or unreadable for %s\n", dmPath.c_str());
        bmClose(fp);
        return 8;
    }
    std::vector<std::pair<uint64_t, uint64_t>> blocks;
    blocks.reserve(fp->idx->nItems);
    for(uint32_t tid = 0; tid < fp->cl->nKeys; ++tid) {
        bmOverlapBlock_t *overlaps = walkRTreeNodes(fp, fp->idx->root, tid, 0, fp->cl->len[tid]);
        if(!overlaps) {
            fprintf(stderr, "[dmcheck] failed to walk index for contig %u in %s\n", tid, dmPath.c_str());
            bmClose(fp);
            return 9;
        }
        for(uint64_t i = 0; i < overlaps->n; ++i) {
            blocks.emplace_back(overlaps->offset[i], overlaps->size[i]);
        }
        destroyOverlapBlock(overlaps);
    }
    if(fp->idx->nItems == 0 || blocks.empty()) {
        fprintf(stderr, "[dmcheck] no indexed data blocks found in %s (records=%" PRIu64 ")\n", dmPath.c_str(), gDmRecordsWritten);
        bmClose(fp);
        return 9;
    }
    for(size_t i = 0; i < blocks.size(); ++i) {
        const uint64_t start = blocks[i].first;
        const uint64_t span = blocks[i].second;
        if(span == 0 || start > fileSize || start + span > fileSize) {
            fprintf(stderr, "[dmcheck] index entry %zu out of bounds in %s (offset=%" PRIu64 ", size=%" PRIu64 ")\n", i, dmPath.c_str(), start, span);
            bmClose(fp);
            return 10;
        }
        if(i && start < blocks[i-1].first) {
            fprintf(stderr, "[dmcheck] index offsets not monotonic at entry %zu in %s\n", i, dmPath.c_str());
            bmClose(fp);
            return 11;
        }
    }
    if(!blocks.empty()) {
        uint8_t headerBuf[32];
        if(bmSetPos(fp, blocks.front().first) || bmRead(headerBuf, 1, sizeof(headerBuf), fp) != sizeof(headerBuf)) {
            fprintf(stderr, "[dmcheck] failed to read first data block header in %s\n", dmPath.c_str());
            bmClose(fp);
            return 12;
        }
        const uint32_t tid = ((uint32_t*)headerBuf)[0];
        const uint32_t start = ((uint32_t*)headerBuf)[1];
        const uint32_t end = ((uint32_t*)headerBuf)[2];
        const uint16_t nItems = ((uint16_t*)headerBuf)[11];
        if(tid >= fp->cl->nKeys || start > end || nItems == 0) {
            fprintf(stderr, "[dmcheck] invalid first data block metadata in %s (tid=%u start=%u end=%u nItems=%u)\n", dmPath.c_str(), tid, start, end, nItems);
            bmClose(fp);
            return 13;
        }
    }

    uint64_t totalIntervals = 0;
    for(uint32_t tid = 0; tid < fp->cl->nKeys; ++tid) {
        const uint32_t chromLen = fp->cl->len[tid];
        const uint32_t window = chromLen > 10000 ? 10000 : chromLen;
        bmOverlappingIntervals_t *hit = bmGetOverlappingIntervals(fp, fp->cl->chrom[tid], 0, window);
        if(hit) {
            totalIntervals += hit->l;
            bmDestroyOverlappingIntervals(hit);
        }
        if(totalIntervals > 0) break;
    }
    if(totalIntervals == 0) {
        fprintf(stderr, "[dmcheck] readable index but no intervals returned from %s (chroms=%u blocks=%zu)\n", dmPath.c_str(), fp->cl->nKeys, blocks.size());
        bmClose(fp);
        return 14;
    }

    bmClose(fp);
    if(verbose) fprintf(stderr, "[dmcheck] %s: ok (%zu indexed blocks, firstIntervals=%" PRIu64 ")\n", dmPath.c_str(), blocks.size(), totalIntervals);
    return 0;
}
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	
	const char* Help_String="Command Format :  dmtools bam2dm [options] -g genome.fa -i/-b <SamfileSorted/BamfileSorted> -m <methratio dm outfile>\n"
		"\nUsage: dmtools bam2dm -g genome.fa -b align.sort.bam -m meth.dm\n"
        "\t [bam2dm] mode paramaters, required\n"
		"\t-g|--genome           genome fasta file\n"
		"\t-b|--binput           Bam format file, sorted by chrom.\n"
        "\t-m|--methratio        methratio.dm output file\n"
		"\t [bam2dm] mode paramaters, options\n"
        "\t-n|--Nmismatch        Number of mismatches, default 0.1 percentage of read length. [0-1]\n"
		"\t-Q                    caculate the methratio while read QulityScore >= Q. default:30\n"
		"\t-c|--coverage         >= <INT> coverage. default:4\n"
		"\t--maxcoverage         <= <INT> coverage. default:1000\n"
		"\t-nC		             >= <INT> nCs per region. default:1\n"
		"\t-r|--remove_dup       REMOVE_DUP, default:false\n"
        "\t--ph                  print head bases of seq, default: seq.txt, motif analysis\n"
        "\t--mrtxt               print prefix.methratio.txt file\n"
        "\t--cf                  context filter for print results, C, CG, CHG, CHH, default: C\n"
        "\t--pe_overlap          skip paired end overlap region, 0 or 1, default 1\n"
        "\t-p|--threads          [int] threads (default: 1)\n"
        "\t--NoMe                data type for NoMe-seq\n"
        "\t [DM format] paramaters\n"
        "\t-C                    print coverage\n"
        "\t-S                    print strand\n"
        "\t--Cx                  print context\n"
        "\t-E                    print end\n"
        "\t--Id                  print ID\n"
        "\t--zl                  The maximum number of zoom levels. [1-10], default: 2\n"
        "\t--chunk-by <chrom|bin>  process by chromosome (default) or fixed-size bins\n"
        "\t--bin-size <INT>      bin size in bases when using --chunk-by bin (default: 2000)\n"
        "\t--debug               verbose logging of parsed options and progress\n"
        "\t--check <dm>          validate an existing dm file and exit\n"
        "\t--validate-output     validate dm after writing (structural check)\n"
        "\t-i|--input            Sam format file, sorted by chrom.\n"
        "\t--countreadC          count number of mC and C for per read\n"
        "\t--chrom               chr/chr:s-e only process this region\n"
        "\t--bedfile             bedfile for process countreadC\n"
        "\t--pms                 file name for print methstate\n"
//        "\t--rrbs                RRBS mode\n"
        "\t-h|--help";

	Char2Comp['A']=Char2Comp['a']='T';
	Char2Comp['C']=Char2Comp['c']='G';
	Char2Comp['G']=Char2Comp['g']='C';
	Char2Comp['T']=Char2Comp['t']='A';
	Char2Comp['N']=Char2Comp['n']='N';
    Char2Comp['U']=Char2Comp['u']='U';
	int Genome_CountX=0;
	char Output_Name[100];
	strcpy(Output_Name, "None");
    char Output_Name_MS[100];
    strcpy(Output_Name_MS, "None");

	string Prefix2="None";
	string methOutfileName;
    string GCHOutfileName;
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
	int zoomlevel = 0;
    char processChr[1000];
    strcpy(processChr, "NAN-mm");
    char processRegion[1000];
    strcpy(processRegion, "NAN-mm");
    char bedfilename[1000];
    strcpy(bedfilename, "");
    std::string binTempBed;
    size_t binTaskCount = 0;
    bool binBedGenerated = false;
    int NTHREADS = 1;
    bool checkOnly = false;
    bool validateOutput = false;
    std::string checkDmFile;
    bool debugMode = false;
    std::string chunkBy = "chrom";
    int binSize = 2000;

	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-f") ||!strcmp(argv[i], "--sam")  )
		{
			strcpy(Output_Name, argv[++i]);
			Sam=1;
		}else if(!strcmp(argv[i], "--pms")  )
        {
            strcpy(Output_Name_MS, argv[++i]);
            printmethstate=1;
        }else if(!strcmp(argv[i], "--rrbs") )
        {
            //rrbs=1;
        }
		else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "--chrom"))
        {
            strcpy(processRegion, argv[++i]);
            char tempRegion[1000];
            strcpy(tempRegion, processRegion);
            strcpy(processChr, strtok(tempRegion, ",:-"));
            fprintf(stderr, "%s %s\n", processRegion, processChr);
        }else if(!strcmp(argv[i], "--bedfile"))
        {
            strcpy(bedfilename, argv[++i]);
        }
        else if(!strcmp(argv[i], "-Q"))
                {
                        QualCut=atoi(argv[++i]);
        }else if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--threads"))
        {
            NTHREADS=atoi(argv[++i]);
            if(NTHREADS < 1) NTHREADS = 1;
        }
        else if(!strcmp(argv[i], "--debug"))
        {
            debugMode = true;
        }
        else if(!strcmp(argv[i], "--chunk-by"))
        {
            if(i + 1 >= argc) {
                fprintf(stderr, "--chunk-by requires chrom or bin\n");
                return 1;
            }
            chunkBy = argv[++i];
            if(chunkBy != "chrom" && chunkBy != "bin") {
                fprintf(stderr, "--chunk-by must be chrom or bin\n");
                return 1;
            }
        }
        else if(!strcmp(argv[i], "--bin-size"))
        {
            if(i + 1 >= argc) {
                fprintf(stderr, "--bin-size requires an integer\n");
                return 1;
            }
            binSize = atoi(argv[++i]);
            if(binSize <= 0) {
                fprintf(stderr, "--bin-size must be positive\n");
                return 1;
            }
        }
        else if(!strcmp(argv[i],"--cf")){
            contextfilter = argv[++i];
        }else if(!strcmp(argv[i],"--NoMe")){
            tech = "NoMe";
        }
        else if(!strcmp(argv[i],"--pe_overlap")){
            skipOverlap = atoi(argv[++i]);
        }else if(!strcmp(argv[i],"--ph")){
            PHead = 1;
            hsPrefix=argv[++i];
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
		}else if(!strcmp(argv[i],"--countreadC"))
        {
            countreadC=1;
        }else if(!strcmp(argv[i],"-o"))
        {
            Methratio=true;
            Prefix2=argv[++i];
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
			//while(i!=(argc-1) && argv[i][0]!='-')
			//{
			//	i++;
			//	continue;
			//}
			//if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
            InFileEnd=InFileStart;
			bamformat=true;
                }
                else if(!strcmp(argv[i], "--check"))
                {
                        checkOnly = true;
			if(i + 1 >= argc) {
				fprintf(stderr, "--check requires a dm filepath\n");
				exit(1);
			}
			checkDmFile = argv[++i];
		}
                else if(!strcmp(argv[i], "--validate-output"))
                {
                        validateOutput = true;
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
    if(DEBUG>1) fprintf(stderr, "\nXXX111\n");
        for(int i = 0; i < argc; i++) {CMD.append(argv[i]); CMD.append(" ");}
    if(strcmp(Prefix.c_str(),"None")==0) Prefix = Prefix2;
    if(!Methratio) onlyPHead = 1;

    if(checkOnly) {
        if(checkDmFile.empty()) {
            fprintf(stderr, "[dmcheck] no dm file provided to --check\n");
            return 1;
        }
        return validateDmFile(checkDmFile, true);
    }
    if(validateOutput && !Methratio) {
        fprintf(stderr, "[dmcheck] --validate-output requires --methratio output path\n");
        return 1;
    }

    gDebugMode = debugMode;
    gDmRecordsWritten = 0;
    gDmRecordsWrittenGch = 0;
    fprintf(stderr, "[bam2dm] Methratio=%d Prefix=%s printtxt=%d chunkBy=%s processChr=%s processRegion=%s\n",
            Methratio ? 1 : 0,
            Prefix.empty() ? "(none)" : Prefix.c_str(),
            printtxt,
            chunkBy.c_str(),
            processChr,
            processRegion);

    if(chunkBy == "bin") {
        if(!bamformat) {
            fprintf(stderr, "--chunk-by bin currently requires BAM input via -b/--binput\n");
            return 1;
        }
        if(!Methratio) {
            fprintf(stderr, "[bin] --chunk-by bin requires methratio output (-m)\n");
            return 1;
        }
        if(InFileStart == 0 || InFileEnd < InFileStart) {
            fprintf(stderr, "[bin] no input BAM provided\n");
            return 1;
        }
        std::string targetChrom = processChr;
        int rc = materializeBinTasksToBed(Geno, targetChrom, binSize, binTempBed, binTaskCount, debugMode);
        if(rc != 0) {
            return rc;
        }
        strncpy(bedfilename, binTempBed.c_str(), sizeof(bedfilename) - 1);
        bedfilename[sizeof(bedfilename) - 1] = '\0';
        binBedGenerated = true;
        strcpy(processChr, "NAN-mm");
    }

    if(debugMode) {
        fprintf(stderr, "[bam2dm] debug on\n");
        fprintf(stderr, "[bam2dm] threads=%d chunk-by=%s bin-size=%d input-format=%s\n", NTHREADS, chunkBy.c_str(), binSize, bamformat ? "bam" : "sam");
        fprintf(stderr, "[bam2dm] genome=%s methratio=%s\n", Geno.empty() ? "(none)" : Geno.c_str(), Prefix.c_str());
        if(binBedGenerated) {
            fprintf(stderr, "[bam2dm] bin tasks materialized to %s (%zu tasks)\n", binTempBed.c_str(), binTaskCount);
        }
    }

    if(chunkBy == "bin") {
        std::vector<BinTask> tasks;
        int rc = loadBinTasksFromBed(binTempBed, tasks, debugMode);
        if(rc != 0) {
            if(binBedGenerated && !binTempBed.empty()) unlink(binTempBed.c_str());
            return rc;
        }
        std::string bamPath = argv[InFileStart];
        std::string partDir = Prefix + ".parts";
        std::vector<BinPartLocator> locators;
        rc = runBinTasksToParts(bamPath, tasks, NTHREADS, partDir, debugMode, locators);
        if(rc == 0) {
            rc = mergeBinPartsToDm(Geno, tasks, locators, std::max(1, NTHREADS), partDir, methOutfileName, zoomlevel, debugMode, validateOutput);
        }
        for(int t = 0; t < std::max(1, NTHREADS); ++t) {
            char partPath[PATH_MAX];
            snprintf(partPath, sizeof(partPath), "%s/thread_%d.tmp", partDir.c_str(), t);
            unlink(partPath);
        }
        rmdir(partDir.c_str());
        if(binBedGenerated && !binTempBed.empty()) unlink(binTempBed.c_str());
        return rc;
    }

        if(argc<=3) printf("%s \n",Help_String);
        if (argc >3  && InFileStart)
        {
		string log;
                log=Prefix;
                log+=".methlog.txt";
		string mCdensity;
                mCdensity=Prefix;  
                mCdensity+=".mCdensity.txt";

                string NCcoverage;
                NCcoverage=Prefix;
                NCcoverage+=".NCcoverage.txt";

		string mCcatero;
                mCcatero=Prefix;
                mCcatero+=".mCcatero.txt";

		fprintf(stderr, "[DM::calmeth] Coverage and validC: %d %d, %d\n", Mcoverage, maxcoverage, nCs);	
		try
		{
			time(&Start_Time);
			Build_Pow10();

			
                        char **chroms = (char **)calloc(MAX_CHROM, sizeof(char*));
                        //if(!chroms) goto error;
                        uint32_t *chrLens = (uint32_t *)malloc(sizeof(uint32_t) * MAX_CHROM);
			
			FILE* GenomeFILE=File_Open(Geno.c_str(),"r");
			fprintf(stderr, "[DM::calmeth] Loading genome sequence : %s\n", Geno.c_str());
                        ARGS args{};
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
						if(chrlen>0) fprintf(stderr, "[DM::calmeth] loaded and store %s %d %d\n", chrName, chrlen, Genome_Count);
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
                
                if(chrlen>0) fprintf(stderr, "[DM::calmeth] loaded and store %s %d %d\n", chrName, chrlen, Genome_Count);

				strcpy(args.Genome_Offsets[Genome_Count].Genome, chrName);
				args.Genome_Offsets[Genome_Count].Offset = chrlen;
				
				chrLens[Genome_Count] = chrlen;
				chroms[Genome_Count] = strdup(chrName);
                chrlen = 0;
                Genome_Count++;
                if(Genome_Count>MAX_CHROM){
                    fprintf(stderr, "\nchrom count bigger than 1000000, please -x-\n");
                    exit(0);
                }
			}

            if(longestChr == 0){
                fprintf(stderr, "Warining: undetected %s in genome file", processChr);
                exit(0);
            }

                        fprintf(stderr, "[DM::calmeth] len count GCcount %lld %d %ld\n", Genome_Size, Genome_Count, totalC+totalG);
            gGenomeSizeBases = static_cast<uint64_t>(Genome_Size);
            if(gGenomeSizeBases > 0) {
                double cgFrac = static_cast<double>(totalC + totalG) / static_cast<double>(gGenomeSizeBases);
                if(cgFrac < 0.01) {
                    fprintf(stderr, "[bam2dm] genome appears converted (C+G fraction %.6f). Use an unconverted reference FASTA.\n", cgFrac);
                    return 1;
                } else if(gDebugMode) {
                    fprintf(stderr, "[bam2dm] genome C+G fraction=%.6f\n", cgFrac);
                }
            }
			//if(!fread(args.Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");
                        if(REMOVE_DUP){
                                args.Marked_Genome=new char[Genome_Size+1];if(!args.Marked_Genome) throw("Insufficient memory to Mark genome..\n");
                                args.Marked_GenomeE=new char[Genome_Size+1];if(!args.Marked_GenomeE) throw("Insufficient memory to Mark genome..\n");
                memset(args.Marked_Genome, 0, Genome_Size+1);
                memset(args.Marked_GenomeE, 0, Genome_Size+1);
                        }
			
			args.OUTFILE = NULL;
            args.OUTFILEMS = NULL;
			assert(longestChr>0);
			fprintf(stderr, "[DM::calmeth] Longest chr: %d\n",longestChr);
			if(Sam && strcmp(Output_Name,"None") ) args.OUTFILE=File_Open(Output_Name,"w");
            if(strcmp(Output_Name_MS,"None") ) args.OUTFILEMS=File_Open(Output_Name_MS,"w");

			char* Split_Point=args.Org_Genome;//split and write...
			args.Genome_List = new Gene_Hash[Genome_Count];
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				String_Hash[args.Genome_Offsets[i].Genome]=i;
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
                                args.Methy_List.plusMethylated = new int[longestChr];
                                args.Methy_List.plusUnMethylated = new int[longestChr];
                                args.Methy_List.NegMethylated = new int[longestChr];
                                args.Methy_List.NegUnMethylated = new int[longestChr];
                                std::fill_n(args.Methy_List.plusG, longestChr, 0);
                                std::fill_n(args.Methy_List.plusA, longestChr, 0);
                                std::fill_n(args.Methy_List.NegG, longestChr, 0);
                                std::fill_n(args.Methy_List.NegA, longestChr, 0);
                                std::fill_n(args.Methy_List.plusMethylated, longestChr, 0);
                                std::fill_n(args.Methy_List.plusUnMethylated, longestChr, 0);
                                std::fill_n(args.Methy_List.NegMethylated, longestChr, 0);
                                std::fill_n(args.Methy_List.NegUnMethylated, longestChr, 0);
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
                                args.methOutFp=File_Open(methOutfileName.c_str(),"w");
                strncpy(args.methOutPath, methOutfileName.c_str(), sizeof(args.methOutPath) - 1);
                args.methOutPath[sizeof(args.methOutPath) - 1] = '\0';
                args.methOutLines = 0;
                if(gDebugMode && args.methOutFp) setvbuf(args.methOutFp, NULL, _IONBF, 0);
                                fprintf(args.methOutFp,"#chromsome\tloci\tstrand\tcontext\tC_count\tCT_count\tmethRatio\teff_CT_count\trev_G_count\trev_GA_count\tMethContext\t5context\n");
                if(gDebugMode) fprintf(stderr, "[dm-writer] mrtxt path: %s fp=%p\n", methOutfileName.c_str(), (void*)args.methOutFp);
                        }

                        fp = NULL;
                        gDmWritersClosed = false;
                        methOutfileName=Prefix;
//                methOutfileName+=".methratio.dm";
                        if(gDebugMode) fprintf(stderr, "[dm-writer] dm output path: %s\n", methOutfileName.c_str());

                                if(bmInit(1<<17) != 0) {
                                        fprintf(stderr, "Received an error in dmInit\n");
                                        return 1;
                                }
                                fp = (binaMethFile_t*)bmOpen((char*)methOutfileName.data(), NULL, "w");
                                if(!fp) {
                                        fprintf(stderr, "An error occurred while opening %s for writing\n", methOutfileName.c_str());
                                        return 1;
                                }
                                fp->type = write_type;

                if(tech=="NoMe"){
                    fp_gch = NULL;
                    GCHOutfileName=Prefix+".gch.dm";

                    if(bmInit(1<<17) != 0) {
                        fprintf(stderr, "Received an error in dmInit for GCH output\n");
                        return 1;
                    }
                    fp_gch = (binaMethFile_t*)bmOpen((char*)GCHOutfileName.data(), NULL, "w");
                    if(!fp_gch) {
                        fprintf(stderr, "An error occurred while opening %s for writing\n", GCHOutfileName.c_str());
                        return 1;
                    }
                    fp_gch->type = write_type;
                }

			}

            args.INbamfilename = new char[1000];
            if(countreadC && strcmp(bedfilename, "")!=0) args.bedFILE=File_Open(bedfilename,"r");

			for(int f=InFileStart;f<=InFileEnd;f++)
			{
				fprintf(stderr, "[DM::calmeth] Processing %d out of %d. File: %s, %d\n\n", f-InFileStart+1,InFileEnd-InFileStart+1, argv[f], bamformat);
				//fseek(args.INFILE, 0L, SEEK_END);args.File_Size=ftello64(args.INFILE);rewind(args.INFILE);
				char s2t[BATBUF];
				//if(args.OUTFILE!=NULL)
 				{ // && printheader
					//samfile_t *bamin = 0;
					args.BamInFile = 0;
					args.header;
                    if(DEBUG>1) fprintf(stderr, "\nXXX3333.1\n");

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
                if(DEBUG>1) fprintf(stderr, "\nXXX3333.2\n");
                                if(f==InFileStart && Methratio){
                                        //Allow up to 10 zoom levels, though fewer will be used in practice
                                        if(bmCreateHdr(fp, zoomlevel)) {
                                                fprintf(stderr, "Failed to create dm header for %s\n", methOutfileName.c_str());
                                                exit(1);
                                        }
                                        //Create the chromosome lists
                                        fp->cl = bmCreateChromList(chroms, chrLens, Genome_Count); //2
                                        if(!fp->cl) {
                                                fprintf(stderr, "Failed to create dm chrom list for %s\n", methOutfileName.c_str());
                                                exit(1);
                                        }
                                        //Write the header
                                        if(bmWriteHdr(fp)) {
                                                fprintf(stderr, "Failed to write dm header for %s\n", methOutfileName.c_str());
                                                exit(1);
                                        }
                                        //Some example methlevel
                                        if(DEBUG>1) fprintf(stderr, "====HHH type %d\n", fp->type);

                        if(tech=="NoMe"){
                        //Allow up to 10 zoom levels, though fewer will be used in practice
                        if(bmCreateHdr(fp_gch, zoomlevel)) {
                                fprintf(stderr, "Failed to create dm header for %s\n", GCHOutfileName.c_str());
                                exit(1);
                        }
                        //Create the chromosome lists
                        fp_gch->cl = bmCreateChromList(chroms, chrLens, Genome_Count); //2
                        if(!fp_gch->cl) {
                                fprintf(stderr, "Failed to create dm chrom list for %s\n", GCHOutfileName.c_str());
                                exit(1);
                        }
                        //Write the header
                        if(bmWriteHdr(fp_gch)) {
                                fprintf(stderr, "Failed to write dm header for %s\n", GCHOutfileName.c_str());
                                exit(1);
                        }
                    }
				}
                if(DEBUG>1) fprintf(stderr, "\nXXX3333.4\n");
				//nothreads
                args.processChr=new char[1000];
                strcpy(args.processChr, processRegion);
                if(DEBUG>1) fprintf(stderr, "\nss000111\n");
				Process_read(&args);
				if(!countreadC) Done_Progress();
				if(!bamformat) fclose(args.samINFILE);
            	if(bamformat)
            	{
                	bam_hdr_destroy(args.header);
                	sam_close(args.BamInFile);
            	}
			}

            if(countreadC && strcmp(bedfilename, "")!=0) fclose(args.bedFILE);
            if(binBedGenerated && !binTempBed.empty()) unlink(binTempBed.c_str());
            if(Sam && strcmp(Output_Name,"None") ) fclose(args.OUTFILE);
            if(printmethstate == 1) fclose(args.OUTFILEMS);

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
                                if(printtxt == 1 && args.methOutFp){
                    if(gDebugMode) fprintf(stderr, "[dm-writer] closing mrtxt %s fp=%p lines=%" PRIu64 "\n", args.methOutPath, (void*)args.methOutFp, args.methOutLines);
                                        fclose(args.methOutFp);
                    args.methOutFp = NULL;
                                }
                         }
			fprintf(stderr, "genome process done!\n");
			
			fprintf(stderr, "Raw count of Met_C in CG:\t%lu\n",met_CG);
			fprintf(stderr, "Raw count of Non_Met_C in CG:\t%lu\n",non_met_CG);
			fprintf(stderr, "Raw count of Met_C in CHG:\t%lu\n",met_CHG);
			fprintf(stderr, "Raw count of Non_Met_C in CHG:\t%lu\n",non_met_CHG);
			fprintf(stderr, "Raw count of Met_C in CHH:\t%lu\n",met_CHH);
			fprintf(stderr, "Raw count of Non_Met_C in CHH:\t%lu\n",non_met_CHH);

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
    myMap.clear();
    bool dmEmpty = Methratio && (gDmRecordsWritten == 0);
    fprintf(stderr, "[DM::calmeth] dm closing\n");
    if(Methratio) closeDmOutputs();
    fprintf(stderr, "[DM::calmeth] dm closed\n");
    if(Methratio && gDebugMode) {
        struct stat dmStat{};
        if(stat(methOutfileName.c_str(), &dmStat) == 0) {
            fprintf(stderr, "[dm-writer] %s size=%lld bytes records=%" PRIu64 "\n", methOutfileName.c_str(),
                    static_cast<long long>(dmStat.st_size), gDmRecordsWritten);
        }
        if(tech=="NoMe") {
            struct stat dmStatGch{};
            if(stat(GCHOutfileName.c_str(), &dmStatGch) == 0) {
                fprintf(stderr, "[dm-writer] %s size=%lld bytes records=%" PRIu64 "\n", GCHOutfileName.c_str(),
                        static_cast<long long>(dmStatGch.st_size), gDmRecordsWrittenGch);
            }
        }
    }
    if(dmEmpty) {
        fprintf(stderr, "[dm-writer] no dm records were written to %s; mark as failure\n", methOutfileName.c_str());
        return 1;
    }
    if(validateOutput && Methratio){
        int validationStatus = validateDmFile(methOutfileName, true);
        if(validationStatus == 0 && tech=="NoMe" && fp_gch){
                    validationStatus = validateDmFile(GCHOutfileName, true);
                }
                if(validationStatus != 0) return validationStatus;
            }
            if(Methratio) bmCleanup();
//delete
                        fprintf(stderr, "[DM::calmeth] Done and release memory!\n");
                        for(int i =0; i < Genome_Count; i++){
                    if(chroms[i]) free(chroms[i]);
            }
                        free(chroms);free(chrLens);
                        if(RELESEM){
                                if(Methratio) freeMethArrays(args);
                                delete [] args.Genome_List; args.Genome_List = NULL;
                                delete [] args.Genome_Offsets; args.Genome_Offsets = NULL;
                                delete [] args.Org_Genome; args.Org_Genome = NULL;
                                delete [] args.Marked_Genome; args.Marked_Genome = NULL;
                        }
		}
		catch(char* Err)
		{
			printf(Err);
			fprintf(stderr, "\nError cigar\n");
			exit(-1);
		}

        fprintf(stderr, "[filter-summary] reads=%" PRIu64 " accepted=%" PRIu64 " bamFiltered=%" PRIu64
                        " mapq=%" PRIu64 " hashMiss=%" PRIu64 " mismatch=%" PRIu64
                        " refOOB=%" PRIu64 " refNonACGT=%" PRIu64 " enterMethyl=%" PRIu64
                        " methylCalls=%" PRIu64 " iterNull=%" PRIu64 "\n",
                gFilterStats.totalReads.load(),
                gFilterStats.processReadAccepted.load(),
                gFilterStats.bamFiltered.load(),
                gFilterStats.mapqFiltered.load(),
                gFilterStats.hashMiss.load(),
                gFilterStats.mismatchFiltered.load(),
                gFilterStats.refOob.load(),
                gFilterStats.refNonACGT.load(),
                gFilterStats.enterMethyl.load(),
                gFilterStats.methylCalls.load(),
                gFilterStats.iterNull.load());

                time(&End_Time);fprintf(stderr, "[DM::calmeth] Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
            exit(0);
            //fclose(OUTLOG);
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

std::string process_cigar(const char* cig,int Read_Len)
{
        char temp[32];
        std::string cigar_rm;
        std::string buffer_cigar;
        unsigned n=0;
        while(*cig!='\0')
        {
                if(*cig>='0' && *cig<='9')
                {
                        if(n+1 >= sizeof(temp)) break;
                        temp[n]=*cig;
                        cig++;n++;
                }else if(*cig=='S')
                {
                        temp[n]='\0';int length=atoi(temp);
                        if(length>0)
                        {
                                buffer_cigar = std::to_string(length) + "S";
                                cigar_rm += buffer_cigar;
                        }
            cig++;n=0;
                }else if(*cig=='M')
                {
                        temp[n]='\0';int length=atoi(temp);
                        if(length>0)
                        {
                                buffer_cigar = std::to_string(length) + "M";
                                cigar_rm += buffer_cigar;
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
                                        cigar_rm = buffer_cigar;
                                }

                        }
                }else if(*cig=='I')
                {
                        temp[n]='\0';int length=atoi(temp);
                        if(length>0)
                        {
                                buffer_cigar = std::to_string(length) + "I";
                                cigar_rm += buffer_cigar;
                        }
                        cig++;n=0;
                }else if(*cig=='D')
                {
                        temp[n]='\0';int length=atoi(temp);
                        if(length>0)
                        {
                                buffer_cigar = std::to_string(length) + "D";
                                cigar_rm += buffer_cigar;
                        }
                        cig++;n=0;
                }else
                {
                        printf(" --%d%c-- ",atoi(temp),*cig);
                continue;//break;
                }
        }
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

//
char **chromsUse_gch;
char **entryid_gch;
uint32_t *starts_gch;
uint32_t *pends_gch;
float *values_gch;
uint16_t *coverages_gch;
uint8_t *strands_gch;
uint8_t *contexts_gch;

static void freeMethArrays(ARGS &args) {
    std::call_once(gMethArraysCleanupFlag, [&](){
        if(args.Methy_List.plusG) { delete[] args.Methy_List.plusG; args.Methy_List.plusG = nullptr; }
        if(args.Methy_List.plusA) { delete[] args.Methy_List.plusA; args.Methy_List.plusA = nullptr; }
        if(args.Methy_List.NegG) { delete[] args.Methy_List.NegG; args.Methy_List.NegG = nullptr; }
        if(args.Methy_List.NegA) { delete[] args.Methy_List.NegA; args.Methy_List.NegA = nullptr; }
        if(args.Methy_List.plusMethylated) { delete[] args.Methy_List.plusMethylated; args.Methy_List.plusMethylated = nullptr; }
        if(args.Methy_List.plusUnMethylated) { delete[] args.Methy_List.plusUnMethylated; args.Methy_List.plusUnMethylated = nullptr; }
        if(args.Methy_List.NegMethylated) { delete[] args.Methy_List.NegMethylated; args.Methy_List.NegMethylated = nullptr; }
        if(args.Methy_List.NegUnMethylated) { delete[] args.Methy_List.NegUnMethylated; args.Methy_List.NegUnMethylated = nullptr; }
    });
}
void print_meth_tofile(int genome_id, ARGS* args){
        if(Methratio)
        {
            fprintf(stderr, "[DM::calmeth] Start process chrom %d\n", genome_id);
            if(gDebugMode) {
                fprintf(stderr, "[dm-writer] chrom=%s dm_fp=%p gch_fp=%p mrtxt_fp=%p mrtxt_lines=%" PRIu64 "\n",
                        args->Genome_Offsets[genome_id].Genome, (void*)fp, (void*)fp_gch, (void*)args->methOutFp, args->methOutLines);
            }
            chromsUse = (char **)calloc(MAX_LINE_PRINT, sizeof(char*));
            entryid = (char **)calloc(MAX_LINE_PRINT, sizeof(char*));
            if(!chromsUse || !entryid) {
                fprintf(stderr, "[dm-writer] failed to allocate output buffers\n");
                return;
            }
            //starts = (uint32_t *)malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
            starts = (uint32_t *)calloc(MAX_LINE_PRINT, sizeof(uint32_t));
            pends = (uint32_t *)malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
            values = (float *)malloc(sizeof(float) * MAX_LINE_PRINT);
            coverages = (uint16_t *)malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
            strands = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
            contexts = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
            const bool enableGch = (tech=="NoMe" && fp_gch);
            if(enableGch) {
                //for GCH chromatin accessibility
                chromsUse_gch = (char **)calloc(MAX_LINE_PRINT, sizeof(char*));
                entryid_gch = (char **)calloc(MAX_LINE_PRINT, sizeof(char*));
                starts_gch = (uint32_t *)calloc(MAX_LINE_PRINT, sizeof(uint32_t));
                pends_gch = (uint32_t *)malloc(sizeof(uint32_t) * MAX_LINE_PRINT);
                values_gch = (float *)malloc(sizeof(float) * MAX_LINE_PRINT);
                coverages_gch = (uint16_t *)malloc(sizeof(uint16_t) * MAX_LINE_PRINT);
                strands_gch = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
                contexts_gch = (uint8_t *)malloc(sizeof(uint8_t) * MAX_LINE_PRINT);
            }

            if(!starts || !pends || !values || !coverages || !strands || !contexts) {
                fprintf(stderr, "[dm-writer] failed to allocate output buffers for dm writing\n");
                return;
            }
            if(enableGch && (!chromsUse_gch || !entryid_gch || !starts_gch || !pends_gch || !values_gch || !coverages_gch || !strands_gch || !contexts_gch)) {
                fprintf(stderr, "[dm-writer] failed to allocate output buffers for GCH dm writing\n");
                return;
            }

        int printL = 0, chrprinHdr = 0, printL_gch = 0, chrprinHdr_gch = 0;
        uint64_t chromRecords = 0;
        uint64_t chromRecordsGch = 0;
        fprintf(stderr, "[DM::calmeth] Processing chrom %d %d %d\n", genome_id, totalC, totalG);
        auto flushMethBuffer = [&](int &len) -> bool {
            if(len == 0) return true;
            int response = chrprinHdr
                ? bmAppendIntervals(fp, starts, pends, values, coverages, strands, contexts, entryid, len)
                : bmAddIntervals(fp, chromsUse, starts, pends, values, coverages, strands, contexts, entryid, len);
            if(response) {
                fprintf(stderr, "bm%sIntervals failed while processing %s (code %d)\n",
                        chrprinHdr ? "Append" : "Add", chromsUse[0], response);
                return false;
            }
            chromRecords += static_cast<uint64_t>(len);
            gDmRecordsWritten += static_cast<uint64_t>(len);
            chrprinHdr = 1;
            for(int i = 0; i < len; ++i) {
                if(chromsUse[i]) {
                    free(chromsUse[i]);
                    chromsUse[i] = NULL;
                }
                starts[i] = 0;
            }
            len = 0;
            return true;
        };
        auto flushGchBuffer = [&](int &len) -> bool {
            if(len == 0) return true;
            if(!enableGch) return true;
            int response = chrprinHdr_gch
                ? bmAppendIntervals(fp_gch, starts_gch, pends_gch, values_gch, coverages_gch, strands_gch, contexts_gch, entryid_gch, len)
                : bmAddIntervals(fp_gch, chromsUse_gch, starts_gch, pends_gch, values_gch, coverages_gch, strands_gch, contexts_gch, entryid_gch, len);
            if(response) {
                fprintf(stderr, "bm%sIntervals failed while processing %s GCH output (code %d)\n",
                        chrprinHdr_gch ? "Append" : "Add", chromsUse_gch[0], response);
                return false;
            }
            chromRecordsGch += static_cast<uint64_t>(len);
            gDmRecordsWrittenGch += static_cast<uint64_t>(len);
            chrprinHdr_gch = 1;
            for(int i = 0; i < len; ++i) {
                if(chromsUse_gch[i]) {
                    free(chromsUse_gch[i]);
                    chromsUse_gch[i] = NULL;
                }
                starts_gch[i] = 0;
            }
            len = 0;
            return true;
        };
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
				//-----F
				std::string context;std::string middleThree;
				
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
                    if(Cover<Mcoverage || Cover > maxcoverage)
                        continue;
                    //methratio
                    float PlusMethratio;
                    if(revGA>0)
                        PlusMethratio=std::min(float(C_count)/(float(C_count+T_count)*revGA),(float)1.0);
                    else
                        PlusMethratio=float(C_count)/float(C_count+T_count);

                    if(tech=="BS-Seq") {
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
					    
                        if(contextfilter != "C") {
                            if(contextfilter != context) continue;
                        }
					    
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
                        if(printL >= MAX_LINE_PRINT) {
                            if(!flushMethBuffer(printL)) goto error;
                        }
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
                        if(printtxt == 1 && args->methOutFp){
                            if(revGA>0) fprintf(args->methOutFp,"%s\t%d\t+\t%s\t%d\t%d\t%f\t%0.001f\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),PlusMethratio,float(C_count+T_count)*revGA,rev_G,(rev_A+rev_G),category.c_str(),Fivecontext.c_str());
                            else fprintf(args->methOutFp,"%s\t%d\t+\t%s\t%d\t%d\t%f\tnull\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),PlusMethratio,rev_G,(rev_A+rev_G),category.c_str(),Fivecontext.c_str());
                            args->methOutLines++;
                            if(gDebugMode && args->methOutLines % 1000000 == 0) {
                                fprintf(stderr, "[dm-writer] mrtxt lines=%" PRIu64 " last chrom=%s\n", args->methOutLines, args->Genome_Offsets[i].Genome);
                            }
                        }
                    } else if(tech=="NoMe") {
                        // 获取中间三个字符的子串
                        middleThree = Fivecontext.substr(1, 3);
                        if(middleThree=="ACG" || middleThree=="TCG") { //WCG for methylation
                            if(printL >= MAX_LINE_PRINT) {
                                if(!flushMethBuffer(printL)) goto error;
                            }
                            chromsUse[printL] = strdup(args->Genome_Offsets[i].Genome);
                            starts[printL] = l+1;
                            pends[printL] = l+2;
                            coverages[printL] = (C_count+T_count);
                            values[printL] = PlusMethratio;
                            strands[printL] = 0; //0 represent '+'
                            contexts[printL] = 0;
                            printL++;
                        }else if(middleThree=="GCA" || middleThree=="GCT" || middleThree=="GCC") { //GCH for chromatin accessibility
                            if(printL_gch >= MAX_LINE_PRINT) {
                                if(enableGch && !flushGchBuffer(printL_gch)) goto error;
                            }
                            chromsUse_gch[printL_gch] = strdup(args->Genome_Offsets[i].Genome);
                            starts_gch[printL_gch] = l+1;
                            pends_gch[printL_gch] = l+2;
                            coverages_gch[printL_gch] = (C_count+T_count);
                            values_gch[printL_gch] = PlusMethratio;
                            strands_gch[printL_gch] = 0; //0 represent '+'
                            contexts_gch[printL_gch] = 0;
                            printL_gch++;
                        }
                    }

				}

				if(args->Methy_List.NegMethylated[l]+args->Methy_List.NegUnMethylated[l]>=1) //coverage neg
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

                    if(Cover<Mcoverage || Cover > maxcoverage) continue;
                    //methratio
                    float NegMethratio;
                    if(revGA>0)
                        NegMethratio=std::min(float(C_count)/(float(C_count+T_count)*revGA),(float)1.0);
                    else
                        NegMethratio=float(C_count)/(float(C_count+T_count));

                    if(tech=="BS-Seq") {
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

                        if(contextfilter != "C") {
                            if(contextfilter != context) continue;
                        }
					    
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
                        if(printL >= MAX_LINE_PRINT) {
                            if(!flushMethBuffer(printL)) goto error;
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
                        if(printtxt == 1 && args->methOutFp){
                            if(revGA>0) fprintf(args->methOutFp,"%s\t%d\t-\t%s\t%d\t%d\t%f\t%0.001f\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),NegMethratio,float(C_count+T_count)*revGA,rev_G,(rev_G+rev_A),category.c_str(),Fcontext);
                            else fprintf(args->methOutFp,"%s\t%d\t-\t%s\t%d\t%d\t%f\tnull\t%d\t%d\t%s\t%s\n",args->Genome_Offsets[i].Genome,l+1,context.c_str(),C_count,(C_count+T_count),NegMethratio,rev_G,(rev_G+rev_A),category.c_str(),Fcontext);
                            args->methOutLines++;
                            if(gDebugMode && args->methOutLines % 1000000 == 0) {
                                fprintf(stderr, "[dm-writer] mrtxt lines=%" PRIu64 " last chrom=%s\n", args->methOutLines, args->Genome_Offsets[i].Genome);
                            }
                        }
                    }else if(tech=="NoMe") {
                        // 获取中间三个字符的子串
                        string fivestring = Fcontext;
                        middleThree = fivestring.substr(1, 3);
                        if(middleThree=="ACG" || middleThree=="TCG") { //WCG for methylation
                            if(printL >= MAX_LINE_PRINT) {
                                if(!flushMethBuffer(printL)) goto error;
                            }
                            chromsUse[printL] = strdup(args->Genome_Offsets[i].Genome);
                            starts[printL] = l+1;
                            pends[printL] = l+2;
                            coverages[printL] = (C_count+T_count);
                            values[printL] = NegMethratio;
                            strands[printL] = 1;
                            contexts[printL] = 0;
                            printL++;
                        }else if(middleThree=="GCA" || middleThree=="GCT" || middleThree=="GCC") { //GCH for chromatin accessibility
                            if(printL_gch >= MAX_LINE_PRINT) {
                                if(enableGch && !flushGchBuffer(printL_gch)) goto error;
                            }
                            chromsUse_gch[printL_gch] = strdup(args->Genome_Offsets[i].Genome);
                            starts_gch[printL_gch] = l+1;
                            pends_gch[printL_gch] = l+2;
                            coverages_gch[printL_gch] = (C_count+T_count);
                            values_gch[printL_gch] = NegMethratio;
                            strands_gch[printL_gch] = 1;
                            contexts_gch[printL_gch] = 0;
                            printL_gch++;
                        }
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
        if(!flushMethBuffer(printL)) goto error;
        if(enableGch) {
            if(!flushGchBuffer(printL_gch)) goto error;
        }
        if(gDebugMode) {
            fprintf(stderr, "[dm-writer] chrom %s records=%" PRIu64 " gch_records=%" PRIu64 "\n", args->Genome_Offsets[genome_id].Genome,
                    chromRecords, chromRecordsGch);
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
        free(chromsUse);
        free(entryid);

        free(starts);
        free(pends); free(values); free(coverages); free(strands); free(contexts);

        if(enableGch) {
            for(i =0; i < MAX_LINE_PRINT; i++){
                if(starts_gch[i]>0 && chromsUse_gch[i]) free(chromsUse_gch[i]);
            }
            fprintf(stderr, "[DM::calmeth] Free mem in all others2\n\n");
            free(chromsUse_gch);
            free(entryid_gch);
            free(starts_gch);
            free(pends_gch); free(values_gch); free(coverages_gch); free(strands_gch); free(contexts_gch);
        }

	}//end methratio
	return;

error:
        closeDmOutputs();
        fprintf(stderr, "\nEEEEE main Received an error somewhere!\n");
        return;
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
        char temp[32];unsigned lens=0;int Glens=0;int RLens=0;
        unsigned n=0;
        int Nmismatch=0;  //chrLen; //((ARGS *)arg)->Genome_Offsets[Hash_Index+1].Offset
        const char* cigr=CIGr.c_str();
        while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
        {
                if(*cigr>='0' && *cigr<='9')
                {
                        if(n+1 >= sizeof(temp)) return Nmismatch;
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
                                if (pos+g-1 >= chrLen) { gFilterStats.refOob.fetch_add(1); break; }
                                char genome_Char = toupper(Genome_seq[pos+g-1]);
                                if(genome_Char!='A' && genome_Char!='C' && genome_Char!='G' && genome_Char!='T') {
                                        gFilterStats.refNonACGT.fetch_add(1);
                                        continue;
                                }
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
                        temp[n]='\0';unsigned int length=atoi(temp);
                        Glens+=length;
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
    if(Flag & 0x8) return -1; //remove later f9o test
    if(!(Flag & 0x2)) return -1; // only keep pair

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
	if(Flag==4 || (int)pos <= 0 ) return -1;
        char strtemp[256];
        size_t j=0;
        if (c->n_cigar == 0) return -1;
        else {
                for (i = 0; i < c->n_cigar; ++i) {
                        int written = snprintf(strtemp, sizeof(strtemp), "%d%c", bam_get_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP=X"[bam_get_cigar(b)[i]&BAM_CIGAR_MASK]);
                        if(written < 0) return -1;
                        if(j + static_cast<size_t>(written) >= BATBUF) {
                                fprintf(stderr, "[bam2dm] CIGAR string too long for buffer on read %s\n", Dummy);
                                return -1;
                        }
                        memcpy(CIG + j, strtemp, static_cast<size_t>(written));
                        j += static_cast<size_t>(written);
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
    FILE* fPH=File_Open(hsPrefix.c_str(),"w");
	int processed_vali = 0;
	int Progress=0;Number_of_Tags=INITIAL_PROGRESS_READS;
	if(!countreadC) Init_Progress();
	int mismatch=0;int pos=0;int Top_Penalty=0;int mapQuality=0;int Flag=-1;
	string readString="";
	int hitType=0;
    int fileprocess = 0;
	
	string hits[MAX_HITS_ALLOWED];
	char Comp_String[MAXTAG];for (int i=1;i<MAXTAG;Comp_String[i++]=0);
	//start to read batman hit file
        char *s2t = (char*) malloc(1000);
        std::vector<char> read_Methyl_Info;
        std::vector<char> rawReadBeforeBS;
        char temp[32];
    char headseq[4]; unsigned int hs=0, hs_r =3; int printh = 0; char headseq_rc[4]; 
	char Dummy[BATBUF],forReadString[BATBUF],Chrom[CHROMSIZE];
	char Chrom_P[CHROMSIZE];int pos_P=0;int Insert_Size=0;int Qsingle=0; //Paired-end reads
	string CIGr;char CIG[BATBUF];
	char forQuality[BATBUF],rcQuality[BATBUF],Quality[BATBUF];
	int bamr=1;char* samaddress;
    struct timeval now;
    struct timespec outtime;
	bam1_t *b = bam_init1();
    hts_idx_t *idx = NULL;
    hts_itr_t *iter;
    int left_end = -1;
    int readC = 0, readmC = 0, readCG = 0, readmCG = 0, readCHG = 0, readmCHG = 0, readCHH = 0, readmCHH = 0;

    char PerLine[2000];
    char chrom[100]; int start=0, end=0;
    char processregion[200];
    PerLine[0] = '\0';
    processregion[0] = '\0';
  int psta = 0;
  while(strcmp(((ARGS *)arg)->processChr, "NAN-mm") != 0 || ( ((ARGS *)arg)->bedFILE != NULL && fgets(PerLine,2000,((ARGS *)arg)->bedFILE)!=NULL) || psta == 0){
    psta = 1;

    if(PerLine[0] != '\0' && PerLine[0] == '#') continue;
    if(PerLine[0] != '\0') {
      sscanf(PerLine, "%s%d%d", chrom, &start, &end);
      sprintf(processregion, "%s:%d-%d", chrom, start, end);
      fprintf(stderr, "Processing... %s\n", processregion);
    }

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
    else if(processregion[0] != '\0') {
        if ((idx = sam_index_load(((ARGS *)arg)->BamInFile, ((ARGS *)arg)->INbamfilename )) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            exit(0);
        }
        if ((iter = sam_itr_querys(idx, ((ARGS *)arg)->header, processregion)) == 0) {
            fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, processregion);
            exit(0);
        }
    }
    bamr=1;

	while( (!bamformat && (samaddress = fgets(s2t,BATBUF,((ARGS *)arg)->samINFILE))!=NULL) || (bamformat && bamr>0 ))
	{

        left_end = -1;

		if(bamr < -1) {
			fprintf(stderr, "\ntruncated file.\n");
		}
		hitType = 0;
		
                Progress++;
                fileprocess++;
                gFilterStats.totalReads.fetch_add(1);

		if(bamformat) 
		{
            if(strcmp(((ARGS *)arg)->processChr, "NAN-mm") == 0 && processregion[0] == '\0'){
                            //(r = samread(( (ARGS *)arg)->BamInFile, b));
                bamr = sam_read1(( (ARGS *)arg)->BamInFile, ((ARGS *)arg)->header, b);
                if(bamr<=0) { if(bamr < 0) gFilterStats.iterNull.fetch_add(1); break; }
                            //bam_tostring(((ARGS *)arg)->header , b, s2t);
                            int ct = processbamread(((ARGS *)arg)->header, b, Dummy,Flag,Chrom,pos,mapQuality,CIG,Chrom_P,pos_P,Insert_Size,forReadString,forQuality, hitType);
                            if(ct == -1) { gFilterStats.bamFiltered.fetch_add(1); continue; }
            }else{
                bamr = sam_itr_next(((ARGS *)arg)->BamInFile, iter, b);
                if(bamr<=0) { if(bamr < 0) gFilterStats.iterNull.fetch_add(1); else gFilterStats.iterNull.fetch_add(1); break; }
                int ct = processbamread(((ARGS *)arg)->header, b, Dummy,Flag,Chrom,pos,mapQuality,CIG,Chrom_P,pos_P,Insert_Size,forReadString,forQuality, hitType);
                if(ct == -1) { gFilterStats.bamFiltered.fetch_add(1); continue; }
            }
                }
        Total_Reads++;

                if ( fileprocess>=1000000  ) {
                        fprintf_time(stderr, "Processed %d reads.\n", Total_Reads);
                        fileprocess = 0;
                        maybeLogFilterStats();
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

        if(mapQuality < QualCut) { gFilterStats.mapqFiltered.fetch_add(1); continue; }

		if(!bamformat){
 			if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) continue;
			if(strcmp(Chrom_P, "*") == 0) {
				pos_P = 0;
				Insert_Size = 0;
			}
		}

                readString=forReadString;
                const size_t cigar_cap = estimateAlignedColumns(CIGr);
                const size_t read_buffer_cap = std::max(readString.size() + 1024, cigar_cap);
                read_Methyl_Info.assign(read_buffer_cap, '\0');
                rawReadBeforeBS.assign(read_buffer_cap, '\0');
                int Read_Len=readString.length();
                CIGr=CIG;
		//for(;forReadString[Read_Len]!=0 && forReadString[Read_Len]!='\n' && forReadString[Read_Len]!='\r';Read_Len++);
	    	iter = String_Hash.find(processingchr.c_str());
		H = -1;
                if(iter != String_Hash.end()){
                        H = iter->second;
                }else { gFilterStats.hashMiss.fetch_add(1); continue; }

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
				hitType=1;  //99
			else if( !(Flag & 0x10)  && (Flag & 0x80) )
				hitType=2; //163
			else if( (Flag & 0x10)  && (Flag & 0x80) )
				hitType=3; //147
			else if( (Flag & 0x10)  && (Flag & 0x40) )
				hitType=4; //83
			if(hitType==0) continue;
            // 1+3, 2+4
		}
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
                if(Nmismatch > 0.5 + UPPER_MAX_MISMATCH * strlen(readString.c_str())) { gFilterStats.mismatchFiltered.fetch_add(1); continue; }

                gFilterStats.processReadAccepted.fetch_add(1);

        if(skipOverlap == 1) {
            if (Flag & 0x1) {
                if(strcmp(Chrom_P, "=") == 0 && pos > pos_P){
                    std::string key = Dummy;
                    auto it = myMap.find(key);
                    if (it != myMap.end()) {
                        left_end = it->second; //myMap[key];
                        myMap.erase(it);
                        //myMap.erase(key);
                    }
                }
            }
        }
		
		unsigned long G_Skip=((ARGS *)arg)->Genome_List[H].Genome - ((ARGS *)arg)->Org_Genome;
            int Flag_rm=0;
		if(hitType == 1 || hitType == 3 ) Flag_rm=2; else Flag_rm=4;
		char Mark = '0';char MarkE='0';
		if(REMOVE_DUP){
            uint64_t markIdx = pos + G_Skip;
            uint64_t markEndIdx = pos + G_Skip + readString.size();
            if(markIdx < gGenomeSizeBases) {
                Mark=((ARGS *)arg)->Marked_Genome[markIdx];
            } else {
                gFilterStats.refOob.fetch_add(1);
                continue;
            }
            if(markEndIdx < gGenomeSizeBases) {
                MarkE=((ARGS *)arg)->Marked_GenomeE[markEndIdx];
            } else {
                gFilterStats.refOob.fetch_add(1);
                continue;
            }
        }
        if( !REMOVE_DUP || (!Mark || !(Mark & Flag_rm)) || (!MarkE || !(MarkE & Flag_rm)) )
		{
                        int Hash_Index=((ARGS *)arg)->Genome_List[H].Index;//load current genome..
                        std::copy(readString.begin(), readString.end(), rawReadBeforeBS.begin());
                        rawReadBeforeBS[readString.size()] = '\0';
                        unsigned lens=0;int Glens=0;int RLens=0;
                        unsigned n=0;bool CONTINUE=false;
                        const char* cigr=CIGr.c_str();
            if(DEBUG>1) fprintf(fPH, "%s %d %s %s\n", Dummy, hitType, cigr, readString.c_str());
            printh = 0;
                        gFilterStats.enterMethyl.fetch_add(1);
                        while(*cigr!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
                        {
                                if(*cigr>='0' && *cigr<='9')
                                {
                                        if(n+1 >= sizeof(temp)) { CONTINUE = true; break; }
                                        temp[n]=*cigr;
                                        cigr++;n++;
				}else if(*cigr=='S')
				{
                                        int i;temp[n]='\0';int length=atoi(temp);
                                        for(i=RLens;i<RLens+length;i++)
                                        {
                                                if(static_cast<size_t>(i) >= read_Methyl_Info.size()) { CONTINUE = true; break; }
                                                read_Methyl_Info[i] = 'S';
                                        }
                                        lens+=length;
                        RLens+=length;
	                cigr++;n=0;
				}else if(*cigr=='M')
				{
                    hs = 0; hs_r = 3;
					temp[n]='\0';int length=atoi(temp);
					for(int k=lens,r=RLens,g=Glens;k<length+lens;r++,k++,g++)
					{
                        //if(strcmp(Dummy, "A00545:105:HWCWLDSX2:3:1105:1118:36182")== 0) fprintf(stderr, "%d %d %d %c, %d %d\n", lens, length, k, readString[k], k-lens, length-4);
						read_Methyl_Info[r] = '=';
						if (pos+g-1 >= ((ARGS *)arg)->Genome_Offsets[Hash_Index].Offset) break;

                        if(pos+g < left_end) {
                            //fprintf(stderr, "AAAAA -------- WWWWW");
                            continue;
                        }

//                        if(rrbs){
//                            if(hitType==1 || hitType == 4) {if(k+2>=RLens) break;}
//                            if(hitType==2 || hitType == 3) {if(k<2) continue;}
//                        }

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
                            if(hitType==1) {
                              if(k-lens<4){
                                if (readString[k]=='T' && genome_Char=='C') headseq[hs] = 'U';
                                else headseq[hs] = readString[k];
                                hs++;
                              }
                            }else{
                              if(k-lens>=length-4) {
                                if (readString[k]=='T' && genome_Char=='C') headseq[hs_r] = 'U';
                                else headseq[hs_r] = readString[k];
                                hs_r--;
                              }
                            }
                            if(PHead == 1) {
                              if(hs == 4 || hs_r == -1) {
                                if(hitType==1 && printh==0) 
                                    fprintf(fPH, "%s\t%d\t%s\t%d\t%s\t%s\n", Dummy, hitType, Chrom, pos, CIGr.c_str(), headseq);
                                printh = 1;
                                if(onlyPHead == 1) break;
                              }
                            }

                                                        if (readString[k]=='C' && genome_Char=='C')
                                                        {
                                readC++; readmC++;
                                                                read_Methyl_Info[r] = 'M';
                                                                if(Methratio ) { gFilterStats.methylCalls.fetch_add(1); ((ARGS *)arg)->Methy_List.plusMethylated[pos+g-1]++; }
                                if(genome_CharFor1=='G')
                                {
                                        readCG++; readmCG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'Z';
                                        met_CG++;
                                }//Z methylated C in CpG context
                                else if(genome_CharFor1!='G' && genome_CharFor2!='G')
                                {
                                        readCHG++; readmCHG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'H';
                                        met_CHH++;
                                }//H methylated C in CHH context
                                else if(genome_CharFor1!='G' && genome_CharFor2=='G')
                                {
                                        readCHH++; readmCHH++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'X';
                                        met_CHG++;
                                }//X methylated C in CHG context
                                                        }
                                                        else if (readString[k]=='T' && genome_Char=='C')
                                                        {
                                readC++;
                                                                read_Methyl_Info[r] = 'U';
                                                                rawReadBeforeBS[k] = 'C';
                                                                if(Methratio ) { gFilterStats.methylCalls.fetch_add(1); ((ARGS *)arg)->Methy_List.plusUnMethylated[pos+g-1]++; }
                                if(genome_CharFor1=='G')
                                {
                                        readCG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'z';
                                        non_met_CG++;
                                }//z unmethylated C in CpG context
                                else if(genome_CharFor1!='G' && genome_CharFor2!='G')
                                {
                                        readCHG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'h';
                                        non_met_CHH++;
                                }//h unmethylated C in CHH context
                                else if(genome_CharFor1!='G' && genome_CharFor2=='G')
                                {
                                        readCHH++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'x';
                                        non_met_CHG++;
                                }//x unmethylated C in CHG context

                                                        }
							else if (readString[k] != genome_Char) 
							{
								read_Methyl_Info[r] = genome_Char;  //readString[k]; // genome_Char; for hypol
							}
							if(Methratio)
							{
								if(readString[k]=='G') ((ARGS *)arg)->Methy_List.plusG[pos+g-1]++;
								else if(readString[k]=='A') ((ARGS *)arg)->Methy_List.plusA[pos+g-1]++;
							}
						}
						else if (hitType==2 || hitType==4) {
                            if(hitType==2) {
                              if(k-lens<4){
                                if (readString[k]=='A' && genome_Char=='G') headseq[hs] = 'U';
                                else headseq[hs] = readString[k];
                                hs++;
                              }
                            }else{
                              if(k-lens>=length-4) {
                                if (readString[k]=='A' && genome_Char=='G') headseq[hs_r] = 'U';
                                else headseq[hs_r] = readString[k];
                                hs_r--;
                              }
                            }
                            if(PHead == 1) {
                              if(hs == 4 || hs_r == -1) {
                                if(hitType==2 && printh == 0) fprintf(fPH, "%s\t%d\t%s\t%d\t%s\t%s\n", Dummy, hitType, Chrom, pos, CIGr.c_str(), headseq);
                                printh = 1;
                                if(onlyPHead == 1) break;
                              }
                            }

                                                        if (readString[k]=='G' && genome_Char=='G')
                                                        {
                                readC++; readmC++;
                                                                read_Methyl_Info[r] = 'M';
                                                                if(Methratio ) { gFilterStats.methylCalls.fetch_add(1); ((ARGS *)arg)->Methy_List.NegMethylated[pos+g-1]++; }
                                if(genome_CharBac1=='C')
                                {
                                        readCG++; readmCG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'Z';
                                        met_CG++;
                                }
                                else if(genome_CharBac1!='C' && genome_CharBac2!='C')
                                {
                                        readCHG++; readmCHG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'H';
                                        met_CHH++;
                                }
                                else if(genome_CharBac1!='C' && genome_CharBac2=='C')
                                {
                                        readCHH++; readmCHH++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'X';
                                        met_CHG++;
                                }
                               }
                                                        else if (readString[k]=='A' && genome_Char=='G')
                                                        {
                                readC++;
                                                                read_Methyl_Info[r] = 'U';
                                                                rawReadBeforeBS[k] = 'G';
                                                                if(Methratio) { gFilterStats.methylCalls.fetch_add(1); ((ARGS *)arg)->Methy_List.NegUnMethylated[pos+g-1]++; }
                                if(genome_CharBac1=='C')
                                {
                                        readCG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'z';
                                        non_met_CG++;
                                }
                                else if(genome_CharBac1!='C' && genome_CharBac2!='C')
                                {
                                        readCHG++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'h';
                                        non_met_CHH++;
                                }
                                else if(genome_CharBac1!='C' && genome_CharBac2=='C')
                                {
                                        readCHH++;
                                                                                if(onlyM == 1) read_Methyl_Info[r] = 'x';
                                        non_met_CHG++;
                                }
                                }
							else if (readString[k] != genome_Char) {
								
								read_Methyl_Info[r] = genome_Char;  //readString[k]; //genome_Char;
							}
							if(Methratio)
							{
								if(readString[k]=='C') ((ARGS *)arg)->Methy_List.NegG[pos+g-1]++;
								else if(readString[k]=='T') ((ARGS *)arg)->Methy_List.NegA[pos+g-1]++;
							}
						}
					}
					cigr++;n=0;lens+=length;Glens+=length;RLens+=length;
                                }else if(*cigr=='I')
                                {
                                        int i;temp[n]='\0';int length=atoi(temp);
                                        for(i=0;i<length;i++)
                                        {
                                                if(static_cast<size_t>(i+RLens) >= read_Methyl_Info.size()) { CONTINUE = true; break; }
                                                read_Methyl_Info[i+RLens] = 'I';
                                        }
                                        lens+=length;
                                        RLens+=length;
                                        cigr++;n=0;
                                }else if(*cigr=='D')
                                {
                                        temp[n]='\0';unsigned int length=atoi(temp);
                                        Glens+=length;
                                        cigr++;n=0;
                                }else
				{
					CONTINUE=true;
					break;
				}
			}
            if((hitType==3 || hitType==4) && printh == 1) {
                onlyComp(headseq_rc, headseq);
                //fprintf(fPH, "%s %d %s %d %s %s\n", Dummy, hitType, Chrom, pos, CIGr.c_str(), headseq_rc);
                fprintf(fPH, "%s\t%d\t%s\t%d\t%s\t%s\n", Dummy, hitType, Chrom, pos, CIGr.c_str(), headseq_rc);
            }
			if(CONTINUE) continue;

            if(skipOverlap == 1) {
                if (Flag & 0x1) {
                    if(strcmp(Chrom_P, "=") == 0 && pos <= pos_P && pos+strlen(readString.c_str()) > pos_P ) {
                        std::string key = Dummy;
                        myMap[key] = pos+Glens;
                    }
                }
            }
			
            if(countreadC){
                if(processregion[0] != '\0')
                 fprintf(stdout ,"%s\t%u\t%s\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",Dummy,Flag,Chrom,pos,readmCG, readCG,readmCHG, readCHG,readmCHH, readCHH, readmC, readC, processregion);
               else if(strcmp(((ARGS *)arg)->processChr, "NAN-mm") != 0)
                 fprintf(stdout ,"%s\t%u\t%s\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",Dummy,Flag,Chrom,pos,readmCG, readCG,readmCHG, readCHG,readmCHH, readCHH, readmC, readC, ((ARGS *)arg)->processChr);
               else fprintf(stdout ,"%s\t%u\t%s\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",Dummy,Flag,Chrom,pos,readmCG, readCG,readmCHG, readCHG,readmCHH, readCHH, readmC, readC);
               readmC=0; readC = 0;readmCG=0; readCG = 0;readmCHG=0; readCHG = 0;readmCHH=0; readCHH = 0;
            }
                        if(!read_Methyl_Info.empty()) {
                                if(static_cast<size_t>(RLens) >= read_Methyl_Info.size()) RLens = static_cast<int>(read_Methyl_Info.size() - 1);
                                read_Methyl_Info[RLens]='\0';
                        }
                        if(!rawReadBeforeBS.empty()) {
                                if(static_cast<size_t>(lens) >= rawReadBeforeBS.size()) lens = static_cast<unsigned>(rawReadBeforeBS.size() - 1);
                                rawReadBeforeBS[lens]='\0';
                        }
                        string mappingStrand="N";
                        if(hitType==1) mappingStrand="YC:Z:CT\tYD:Z:f";
			else if(hitType==4) mappingStrand="YC:Z:CT\tYD:Z:r";
			else if(hitType==3) mappingStrand="YC:Z:GA\tYD:Z:r";
			else if(hitType==2) mappingStrand="YC:Z:GA\tYD:Z:f";
            if(printmethstate == 1 && ((ARGS *)arg)->OUTFILEMS != NULL) {
                flockfile(((ARGS *)arg)->OUTFILEMS);
                fprintf(((ARGS *)arg)->OUTFILEMS, "%s\t%s\t%s\t%s\t%s\t%u\t%d\n", Dummy, read_Methyl_Info.data(), readString.c_str(), CIGr.c_str(), Chrom, pos, hitType);
                funlockfile(((ARGS *)arg)->OUTFILEMS);
            }
            if(Sam && ((ARGS *)arg)->OUTFILE != NULL) {
			    flockfile(((ARGS *)arg)->OUTFILE);
			    if(Sam ) //&& Nmismatch <= UPPER_MAX_MISMATCH )
			    {
				if(SamSeqBeforeBS) 
				{
                                        fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\tRA:Z:%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,rawReadBeforeBS.data(),forQuality,Nmismatch,read_Methyl_Info.data(),mappingStrand.c_str(), readString.c_str());
                                }else
                                {
                                        fprintf(((ARGS *)arg)->OUTFILE,"%s\t%u\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\t%s\n",Dummy,Flag,Chrom,pos,mapQuality,CIGr.c_str(),Chrom_P,pos_P,Insert_Size,readString.c_str(),forQuality,Nmismatch,read_Methyl_Info.data(),mappingStrand.c_str());
                                }
			    }
			    funlockfile(((ARGS *)arg)->OUTFILE);
            }
		}
		if(REMOVE_DUP){
                   const uint64_t markIdx2 = pos+G_Skip;
                   const uint64_t markEndIdx2 = pos+G_Skip+readString.size();
                   if(markIdx2 < gGenomeSizeBases) ((ARGS *)arg)->Marked_Genome[markIdx2] |= Flag_rm; else gFilterStats.refOob.fetch_add(1);
                   if(markEndIdx2 < gGenomeSizeBases) ((ARGS *)arg)->Marked_GenomeE[markEndIdx2] |= Flag_rm; else gFilterStats.refOob.fetch_add(1);
		}
	}

    if(((ARGS *)arg)->bedFILE == NULL || strcmp(((ARGS *)arg)->processChr, "NAN-mm") != 0 ) break;
  }
        if(gDebugMode) {
                maybeLogFilterStats();
        }
        free(s2t);
	fclose(fIS);
    fclose(fPH);
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

float Pr(float Q)
{
        assert((int)Q>=0 && (int)Q<POWLIMIT-1);
        return(1-POW10[(int)Q]);
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

string remove_soft_split(string & cigar,int & Read_Len,int & pos)
{

	const char* cig=cigar.c_str();
        char temp[32];
	int n=0,lens=0,length=0,RLens=0,Glens=0;
	int genome_move_size=0; int headClip=0;int moveSize=0;
	char cigar_rm[1000];*cigar_rm=0;
	char buffer_cigar[100];*buffer_cigar=0; bool cigar_remove=false;
	char upper='N';char upper_cigar[1024]="\0";int upper_length=0;
	while(*cig!='\0')//RLens--READs Length \\ lens--raw reads length \\ GLens--genome Lens
	{
		if(*cig>='0' && *cig<='9')
		{
                        if(n+1 >= sizeof(temp)) break;
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


void Reverse_Comp(char* Dest, char* seq)
{
    int stringlength = strlen(seq);
    if(stringlength<=0) {
        fprintf(stderr, "\nError string %d %d %s\n", stringlength, seq);
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


void onlyComp(char* Dest, char* seq)
{
    int stringlength = strlen(seq);
    if(stringlength<=0) {
        fprintf(stderr, "\nError string %d %d %s\n", stringlength, seq);
        exit(0);
    }

        for (int i=0;i<stringlength;i++)
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
