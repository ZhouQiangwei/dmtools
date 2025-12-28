#include "binaMeth.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int write_dm(const char *path, int use_quant, uint32_t scale, const char **chroms, uint32_t *lens,
                    const char **chromsUse, uint32_t *starts, uint32_t *ends, float *values, uint16_t *coverage,
                    size_t n) {
    if(bmInit(1 << 17) != 0) return 1;
    binaMethFile_t *fp = bmOpen((char *)path, NULL, "w");
    if(!fp) return 2;
    fp->type = BM_MAGIC | BM_END | BM_COVER;
    if(use_quant) fp->type |= BM_VAL_U16;
    fp->valScale = use_quant ? scale : 0;
    fp->valEncoding = use_quant ? 1 : 0;
    if(bmCreateHdr(fp, 0)) { bmClose(fp); return 3; }
    fp->cl = bmCreateChromList((char **)chroms, lens, 1);
    if(!fp->cl) { bmClose(fp); return 4; }
    if(bmWriteHdr(fp)) { bmClose(fp); return 5; }

    if(bmAddIntervals(fp, (char **)chromsUse, starts, ends, values, coverage, NULL, NULL, NULL, (int)n)) {
        bmClose(fp);
        return 6;
    }
    bmClose(fp);
    return 0;
}

static int validate_with_cli(const char *dmPath) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "./dmtools validate -i %s > /dev/null", dmPath);
    return system(cmd);
}

int main(void) {
    const char *chroms[] = {"chrQ"};
    uint32_t lens[] = {1000};
    const char *chromsUse[] = {"chrQ", "chrQ", "chrQ", "chrQ"};
    uint32_t starts[] = {1, 100, 250, 900};
    uint32_t ends[] = {2, 150, 300, 901};
    float values[] = {0.05f, 0.25f, 0.73f, 1.0f};
    uint16_t coverage[] = {5, 10, 20, 30};
    const size_t n = sizeof(values) / sizeof(values[0]);
    const uint32_t scale = 1000;

    const char *float_dm = "test/quant_float.dm";
    const char *quant_dm = "test/quant_u16.dm";

    if(write_dm(float_dm, 0, 0, chroms, lens, chromsUse, starts, ends, values, coverage, n) != 0) return 10;
    if(write_dm(quant_dm, 1, scale, chroms, lens, chromsUse, starts, ends, values, coverage, n) != 0) return 11;

    if(validate_with_cli(float_dm) != 0) return 12;
    if(validate_with_cli(quant_dm) != 0) return 13;

    binaMethFile_t *fp = bmOpen((char *)quant_dm, NULL, "r");
    if(!fp) return 14;
    bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, (char *)chroms[0], 0, lens[0]);
    if(!o || o->l != n) { bmClose(fp); return 15; }

    float tol = 1.0f / (float)scale + 1e-6f;
    for(size_t i = 0; i < n; i++) {
        if(o->start[i] != starts[i] || o->end[i] != ends[i]) { bmDestroyOverlappingIntervals(o); bmClose(fp); return 16; }
        if(o->coverage[i] != coverage[i]) { bmDestroyOverlappingIntervals(o); bmClose(fp); return 17; }
        float err = fabsf(o->value[i] - values[i]);
        if(err > tol) { bmDestroyOverlappingIntervals(o); bmClose(fp); return 18; }
    }
    bmDestroyOverlappingIntervals(o);
    bmClose(fp);
    bmCleanupOnce();
    remove(float_dm);
    remove(quant_dm);
    return 0;
}
