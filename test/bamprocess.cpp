#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include "htslib/sam.h"

enum filetype {
	FBAM = 1,   // BAM file
	FSAM = 2,   // SAM file
};

int ftype;

/*
@discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
*/

string getCigar(const bam1_t *b) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    stringstream ss;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << len << op;
    }
    return ss.str();
}

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}

// get seq 序列的信息记录在8bit 的数据结构中，前4bit 是前面的碱基，后4bit 是后面的碱基
string getSeq(const bam1_t *b) {
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF); 
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}

// get seq quality
string getQual(const bam1_t *b) {
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        s[i] = (char)(data[i] + 33); // 转换成打印的ascci
    }
    return s;
}

// get align aux  
//string md = getAux(aln,"MD");
string getAux(const bam1_t *b,const char tag[2])
{
    kstring_t res = KS_INITIALIZE;    // 需要初始化
    if(bam_aux_get_str(b,tag,&res) == 1) //kstring的string buffer 没有\0终止
    {
        int len = ks_len(&res);
        char *ks_s = ks_str(&res);
        string s(len, '\0');
        for (int i = 0;i<len;i++ ){
            s[i] = ks_s[i];
        }
        ks_free(&res); // 释放资源
        return s;
    }    
    else 
    {
        cerr << "no tag :" << tag << '\n';
        ks_free(&res);
        return "";
    }
    
}

int main_read_per(int argc,char *argv[])
{
    if(argc < 2) {
        cerr << "need bam path";
        return -1;
    }
    samFile *bam_in = sam_open(argv[1],"r"); // open bam file
    bam_hdr_t *bam_header = sam_hdr_read(bam_in); // read header
    bam1_t *aln = NULL;
    aln = bam_init1(); //initialize an alignment
    int hitType = 0;
    int Flag = 0;
    int pos = 0;
    int mapQuality = 0;
    string queryname = "";
    string cigar = "";
    string seq = "";
    string qual = "";

    //return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error
    while (sam_read1(bam_in,bam_header,aln) >= 0) 
    {
        int pos = aln->core.pos ;
        string chr = "*";
        if (aln->core.tid != -1) // 存在无法比对到基因组的reads
            chr = bam_header->target_name[aln->core.tid]; // config name(chromosome)   
        Flag = aln->core.flag;
        if(aln->core.tid == -1) {
            continue;
        }
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

        pos = aln->core.pos + 1;
    	mapQuality = aln->core.qual;

        if(mapQuality < QualCut || Flag==4 || (int)pos <= 0 ) continue;

        // mate chrom
        // c-> replace by aln->core.
        if (aln->core.mtid < 0) Chrom_P[0] =  '*';
        else if (aln->core.mtid == aln->core.tid) sprintf(Chrom_P, "="); 
        else {
            strcpy(Chrom_P, bam_header->target_name[aln->core.mtid];
        }
        //
        if(strcmp(Chrom_P, "*") == 0) {
        	pos_P = 0;
            Insert_Size = 0;
        }else{
		    pos_P = c->mpos + 1;
		    Insert_Size=c->isize;
	    }

        queryname = bam_get_qname(aln);
        cigar = getCigar(aln);
        seq = getSeq(aln);
        qual = getQual(aln);

        cout << "QueryName: " << queryname << '\n' 
            << "Positon: " << chr << '\t' <<  pos <<'\n'
            << "Cigar: " << cigar << '\n'
            << "Seq: " << seq << '\n'
            << "Qual: " << qual << '\n';
        
    }

    bam_destroy1(aln); // 回收资源
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
    cout << "Finshed!\n";
    return 0;
}

int main_region(int argc,char *argv[])
{
    if(argc < 3) {
        cerr << "need bam path; chrom: start pos-end pos";
        return -1;
    }
    samFile *bam_in = sam_open(argv[1],"r"); // open bam file
    hts_idx_t *bam_index = sam_index_load(bam_in,argv[1]); //load index 
    bam_hdr_t *bam_header = sam_hdr_read(bam_in); // read header
    bam1_t *aln = NULL;
    aln = bam_init1(); //initialize an alignment
    /*
    Regions are parsed by hts_parse_reg(), and take one of the following forms:
    region          | Outputs
    --------------- | -------------
    REF             | All reads with RNAME REF
    REF:            | All reads with RNAME REF
    REF:START       | Reads with RNAME REF overlapping START to end of REF
    REF:-END        | Reads with RNAME REF overlapping start of REF to END
    REF:START-END   | Reads with RNAME REF overlapping START to END
    .               | All reads from the start of the file
    *               | Unmapped reads at the end of the file (RNAME '*' in SAM)
    */
        string regin = argv[2]; 
       
    hts_itr_t *iter = sam_itr_querys(bam_index,bam_header,regin);
    if(!iter) cerr << "invalid regin\n";
   // while (sam_read1(bam_in,bam_header,aln) >= 0)
    while (sam_itr_next(bam_in, iter, aln) >= 0)
    {
        int pos = aln->core.pos ;
        string chr = "*";
        if (aln->core.tid != -1)
            chr = bam_header->target_name[aln->core.tid]; // config name(chromosome)            
        string queryname = bam_get_qname(aln);
        string cigar = getCigar(aln);
        string seq = getSeq(aln);
        string qual = getQual(aln);
        string md = getAux(aln,"MD");

        //int64_t newPos = pos - 5;
        cout << "QueryName: " << queryname << '\n' 
            << "Positon: " <<aln->core.tid << '\t' << chr << '\t' <<  pos <<'\n'
            << "Cigar: " << cigar << '\n'
            << md << '\n'
            << "Seq: " << seq << '\n'
            << "Qual: " << qual << '\n';
        
    }
    sam_itr_destroy(iter) ;// free iter
    bam_destroy1(aln);
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
    cout << "Finshed!\n";
    return 0;
}


int sam_test_extract(int argc, char **argv, int optind, htsFile *in, htsFile *out) {
		sam_hdr_t *hdr;
		bam1_t *b;
		hts_idx_t *idx = NULL;
		hts_itr_t *iter = NULL;
		int ret;

		if ((hdr = sam_hdr_read(in)) == NULL) {
			fprintf(stderr, "[E::%s] couldn't read header for '%s'\n", __func__, argv[optind]);
			return  -1;
		}
		if ((b = bam_init1()) == NULL) {
			fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
			goto fail;
		}
		if (ftype == FBAM && optind + 2 <= argc) { // BAM input and has a region.
			if ((idx = sam_index_load(in, argv[optind])) == 0) {
				fprintf(stderr, "[E::%s] fail to load the index for '%s'\n", __func__, argv[optind]);
				goto fail;
			}
			if ((iter = sam_itr_querys(idx, hdr, argv[optind + 1])) == 0) {
				fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[optind + 1]);
				goto fail;
			}
			while ((ret = sam_itr_next(in, iter, b)) >= 0) {
				if (sam_write1(out, hdr, b) < 0) {
					fprintf(stderr, "[E::%s] Error writing output.\n", __func__);
					goto fail;
				}
			}
			if (ret < -1) {
				fprintf(stderr, "[E::%s] Error reading input.\n", __func__);
				goto fail;
			}
			hts_itr_destroy(iter);
			iter = NULL;
			hts_idx_destroy(idx);
			idx = NULL;
		} else if (optind + 2 > argc) {
			while ((ret = sam_read1(in, hdr, b)) >= 0) {
				if (sam_write1(out, hdr, b) < 0) {
					fprintf(stderr, "[E::%s] Error writing alignments.\n", __func__);
					goto fail;
				}
			}
			if (ret < -1) {
				fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
				goto fail;
			}
		} else { // SAM input and has a region.
			fprintf(stderr, "[E::%s] couldn't extract alignments directly from raw sam file.\n", __func__);
			goto fail;
		}
		bam_destroy1(b);
		sam_hdr_destroy(hdr);
		return  0;
	fail:
		if (iter) sam_itr_destroy(iter);
		if (b) bam_destroy1(b);
		if (idx) hts_idx_destroy(idx);
		if (hdr) sam_hdr_destroy(hdr);
		return  1;
}

int main_only(int argc, char **argv) {
	htsFile *in, *out;
	int c, ret, exit_code;
	char moder[8];
	//char modew[800];
	char *outfn = "-";

	ftype = FSAM;
	exit_code = 0;
	strcpy(moder, "r");
	while ((c = getopt(argc, argv, "bo:")) >= 0) {
		switch (c) {
			case 'b': strcat(moder, "b"); ftype = FBAM; break;
			case 'o': outfn = optarg; break;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: %s [-b] [-o out.sam] <in.bam>|<in.sam> [region]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "\t-b:\tUse BAM as input if this option is set, otherwise use SAM as input.\n");
		fprintf(stderr, "\t-o:\tPath to the output file. Output to stdout if this option is not set.\n");
		return  -1;
	}
	if ((in = hts_open(argv[optind], moder)) == NULL) {
		fprintf(stderr, "Error opening '%s'\n", argv[1]);
		return  -3;
	}
	if ((out = hts_open(outfn, "w")) == NULL) {
		fprintf(stderr, "Error opening '%s'\n", argv[2]);
		return  -3;
	}
	if ((ret = sam_test_extract(argc, argv, optind, in, out)) != 0) {
		fprintf(stderr, "Error extracting alignment from '%s'\n", argv[optind]);
		exit_code = -5;
	}
	if ((ret = hts_close(out)) < 0) {
		fprintf(stderr, "Error closing output.\n");
		exit_code = -3;
	}
	if ((ret = hts_close(in)) < 0) {
		fprintf(stderr, "Error closing input.\n");
		exit_code = -3;
	}
	return exit_code;
}