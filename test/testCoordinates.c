#include "binaMeth.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    const char *dmPath = "test/coords.dm";
    if(bmInit(1<<17) != 0) {
        fprintf(stderr, "bmInit failed\n");
        return 1;
    }

    binaMethFile_t *fp = bmOpen((char*)dmPath, NULL, "w");
    if(!fp) {
        fprintf(stderr, "bmOpen write failed\n");
        return 1;
    }
    fp->type = BM_MAGIC | BM_END | BM_COVER | BM_STRAND | BM_CONTEXT;
    if(bmCreateHdr(fp, 0)) {
        fprintf(stderr, "bmCreateHdr failed\n");
        bmClose(fp);
        return 1;
    }

    char *chroms[] = {"chr1"};
    uint32_t lens[] = {10};
    fp->cl = bmCreateChromList(chroms, lens, 1);
    if(!fp->cl) {
        fprintf(stderr, "bmCreateChromList failed\n");
        bmClose(fp);
        return 1;
    }
    if(bmWriteHdr(fp)) {
        fprintf(stderr, "bmWriteHdr failed\n");
        bmClose(fp);
        return 1;
    }

    char *chromsUse[] = {"chr1", "chr1", "chr1"};
    uint32_t starts[] = {1, 2, 5};
    uint32_t ends[] = {2, 3, 6};
    float values[] = {0.1f, 0.2f, 0.3f};
    uint16_t coverage[] = {10, 20, 30};
    uint8_t strands[] = {0, 1, 0};
    uint8_t contexts[] = {1, 2, 3};

    if(bmAddIntervals(fp, chromsUse, starts, ends, values, coverage, strands, contexts, NULL, 3)) {
        fprintf(stderr, "bmAddIntervals failed\n");
        bmClose(fp);
        return 1;
    }
    bmClose(fp);

    fp = bmOpen((char*)dmPath, NULL, "r");
    if(!fp) {
        fprintf(stderr, "bmOpen read failed\n");
        return 1;
    }
    bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, (char*)"chr1", 1, 6);
    if(!o || o->l < 3) {
        fprintf(stderr, "bmGetOverlappingIntervals returned insufficient records\n");
        bmClose(fp);
        return 1;
    }
    if(o->start[0] != 1 || o->end[0] != 2 || o->start[1] != 2 || o->end[1] != 3 || o->start[2] != 5 || o->end[2] != 6) {
        fprintf(stderr, "coordinate mismatch: %u-%u %u-%u %u-%u\n",
                o->start[0], o->end[0], o->start[1], o->end[1], o->start[2], o->end[2]);
        bmDestroyOverlappingIntervals(o);
        bmClose(fp);
        return 1;
    }
    bmDestroyOverlappingIntervals(o);
    bmClose(fp);
    bmCleanupOnce();
    remove(dmPath);
    return 0;
}
