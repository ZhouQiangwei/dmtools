#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"

static int write_fasta(const char *path) {
    FILE *fp = fopen(path, "w");
    if(!fp) return 1;
    fprintf(fp, ">chrA\n");
    fprintf(fp, "CCCCGGGGCC\n");
    fprintf(fp, ">chrB\n");
    fprintf(fp, "TTTTTTTTTT\n");
    fclose(fp);
    return fai_build(path);
}

static int write_bam(const char *path) {
    samFile *fp = sam_open(path, "wb");
    if(!fp) return 1;

    sam_hdr_t *hdr = sam_hdr_init();
    if(!hdr) return 1;
    if(sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "coordinate", NULL) != 0) return 1;
    if(sam_hdr_add_line(hdr, "SQ", "SN", "chrB", "LN", "10", NULL) != 0) return 1;
    if(sam_hdr_add_line(hdr, "SQ", "SN", "chrA", "LN", "10", NULL) != 0) return 1;
    if(sam_hdr_write(fp, hdr) != 0) return 1;

    bam1_t *b = bam_init1();
    if(!b) return 1;

    uint32_t cigar = bam_cigar_gen(4, BAM_CMATCH);
    const char *qual = "IIII";

    if(bam_set1(b, 5, "read1", 99, 1, 0, 60,
                1, &cigar, 1, 4, 0, 4, "CCCC", qual, 0) < 0) return 1;
    if(sam_write1(fp, hdr, b) < 0) return 1;

    if(bam_set1(b, 5, "read2", 147, 1, 4, 60,
                1, &cigar, 1, 0, 0, 4, "GGGG", qual, 0) < 0) return 1;
    if(sam_write1(fp, hdr, b) < 0) return 1;

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    if(sam_close(fp) != 0) return 1;

    if(sam_index_build(path, 0) != 0) return 1;
    return 0;
}

int main(int argc, char **argv) {
    if(argc < 3) {
        fprintf(stderr, "usage: %s <fasta> <bam>\n", argv[0]);
        return 1;
    }
    if(write_fasta(argv[1]) != 0) {
        fprintf(stderr, "failed to write fasta %s\n", argv[1]);
        return 1;
    }
    if(write_bam(argv[2]) != 0) {
        fprintf(stderr, "failed to write bam %s\n", argv[2]);
        return 1;
    }
    return 0;
}
