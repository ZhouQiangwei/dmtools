#include "binaMeth.h"
#include "dmCommon.h"
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <sys/wait.h>
#include <unistd.h>

extern const char* Help_String_scaggregate;

typedef struct {
    char **names;
    size_t n;
} dm_sc_idmap_t;

typedef struct {
    uint64_t n_sites;
    uint64_t total_cov;
    double sum_meth_weighted;
    uint64_t n_sites_cg;
    uint64_t n_sites_chh;
} dm_sc_qc_entry_t;

static void dm_sc_free_idmap(dm_sc_idmap_t *map) {
    if (!map || !map->names) return;
    for (size_t i = 0; i < map->n; i++) free(map->names[i]);
    free(map->names);
    map->names = NULL;
    map->n = 0;
}

static int dm_sc_load_idmap(const char *dmPath, dm_sc_idmap_t *map) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s.idmap.tsv", dmPath);
    FILE *f = fopen(path, "r");
    if (!f) {
        snprintf(path, sizeof(path), "%s.idmap", dmPath);
        f = fopen(path, "r");
    }
    if (!f) {
        fprintf(stderr, "Error: missing idmap for %s (expected %s.idmap.tsv)\n", dmPath, dmPath);
        return 1;
    }
    char line[4096];
    size_t cap = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *idStr = strtok(line, "\t\n");
        char *name = strtok(NULL, "\t\n");
        if (!idStr || !name) continue;
        uint32_t id = (uint32_t)strtoul(idStr, NULL, 10);
        if (id == 0) continue;
        if (id >= map->n) {
            size_t newN = id + 1;
            if (newN > cap) {
                cap = cap ? cap * 2 : 64;
                if (cap < newN) cap = newN;
                char **tmp = realloc(map->names, cap * sizeof(char *));
                if (!tmp) { fclose(f); return 1; }
                for (size_t i = map->n; i < cap; i++) tmp[i] = NULL;
                map->names = tmp;
            }
            map->n = newN;
        }
        if (!map->names[id]) map->names[id] = strdup(name);
    }
    fclose(f);
    if (map->n == 0) return 1;
    return 0;
}

static void dm_sc_qc_free_entries(dm_sc_qc_entry_t *entries) {
    if (!entries) return;
    free(entries);
}

static int dm_sc_qc_parse_context(const char *ctx, uint8_t *val) {
    if (strcasecmp(ctx, "C") == 0) {
        *val = 0;
        return 0;
    } else if (strcasecmp(ctx, "CG") == 0) {
        *val = 1;
        return 0;
    } else if (strcasecmp(ctx, "CHG") == 0) {
        *val = 2;
        return 0;
    } else if (strcasecmp(ctx, "CHH") == 0) {
        *val = 3;
        return 0;
    }
    return -1;
}


static void dm_sc_qc_usage() {
    fprintf(stderr,
            "Usage: dmtools sc-qc -i <input.dm> [-i <input2.dm> ...] -o <output.tsv> [options]\n"
            "Options:\n"
            "  -i, --input <dm>         Input DM file with IDs (can be repeated)\n"
            "  -o, --output <tsv>       Output TSV file for QC metrics\n"
            "      --context <CTX>      Context filter: C, CG, CHG, CHH (default: no filter)\n"
            "      --min-coverage <N>   Minimum coverage per site (default: 1)\n"
            "  -h, --help               Show this help message\n"
            "\nOutput columns:\n"
            "  cell_id  n_sites  total_coverage  mean_coverage  mean_meth  n_sites_cg  n_sites_chh  cg_chh_ratio\n");
}

int dm_sc_qc_main(int argc, char **argv) {
    char **inputs = NULL;
    size_t nInputs = 0, capInputs = 0;
    char *outputPath = NULL;
    uint8_t contextFilter = 0;
    int contextFilterSet = 0;
    uint64_t minCoverage = 1;

    static struct option long_opts[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"context", required_argument, 0, 1},
        {"min-coverage", required_argument, 0, 2},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt, long_idx = 0;
    while ((opt = getopt_long(argc, argv, "hi:o:", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 'i': {
                if (nInputs == capInputs) {
                    size_t newCap = capInputs ? capInputs * 2 : 4;
                    char **tmp = realloc(inputs, newCap * sizeof(char *));
                    if (!tmp) {
                        fprintf(stderr, "Error: memory allocation failed for inputs.\n");
                        goto error;
                    }
                    inputs = tmp;
                    capInputs = newCap;
                }
                inputs[nInputs++] = strdup(optarg);
                if (!inputs[nInputs - 1]) {
                    fprintf(stderr, "Error: memory allocation failed for input path.\n");
                    goto error;
                }
                break;
            }
            case 'o':
                outputPath = strdup(optarg);
                if (!outputPath) {
                    fprintf(stderr, "Error: memory allocation failed for output path.\n");
                    goto error;
                }
                break;
            case 1: /* --context */
                if (dm_sc_qc_parse_context(optarg, &contextFilter) != 0) {
                    fprintf(stderr, "Error: invalid context '%s'. Use C, CG, CHG, or CHH.\n", optarg);
                    goto error;
                }
                contextFilterSet = 1;
                break;
            case 2: /* --min-coverage */
                minCoverage = strtoull(optarg, NULL, 10);
                break;
            case 'h':
                dm_sc_qc_usage();
                goto success;
            default:
                dm_sc_qc_usage();
                goto error;
        }
    }

    if (nInputs == 0 || !outputPath) {
        fprintf(stderr, "Error: -i <input.dm> and -o <output.tsv> are required.\n");
        dm_sc_qc_usage();
        goto error;
    }

    dm_sc_idmap_t idmap = {0};
    dm_sc_qc_entry_t *entries = NULL;

    for (size_t fi = 0; fi < nInputs; fi++) {
        binaMethFile_t *fp = bmOpen(inputs[fi], NULL, "r");
        if (!fp) {
            fprintf(stderr, "Error: cannot open input file %s\n", inputs[fi]);
            goto error_entries;
        }
        if (!(fp->hdr->version & BM_ID)) {
            fprintf(stderr, "Error: sc-qc requires dm files with ID (BM_ID set). File %s has no ID field.\n", inputs[fi]);
            bmClose(fp);
            goto error_entries;
        }
        if (!(fp->hdr->version & BM_COVER)) {
            fprintf(stderr, "Error: sc-qc requires coverage information. File %s lacks BM_COVER.\n", inputs[fi]);
            bmClose(fp);
            goto error_entries;
        }
        if (contextFilterSet && !(fp->hdr->version & BM_CONTEXT)) {
            fprintf(stderr, "Error: context filter requested but file %s lacks BM_CONTEXT.\n", inputs[fi]);
            bmClose(fp);
            goto error_entries;
        }

        if (idmap.n == 0) {
            if (dm_sc_load_idmap(inputs[fi], &idmap) != 0) {
                bmClose(fp);
                goto error_entries;
            }
            entries = calloc(idmap.n, sizeof(dm_sc_qc_entry_t));
            if (!entries) { bmClose(fp); goto error_entries; }
        }

        for (uint64_t ci = 0; ci < fp->cl->nKeys; ci++) {
            char *chrom = fp->cl->chrom[ci];
            uint32_t chromLen = fp->cl->len[ci];
            bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, chrom, 0, chromLen);
            if (!o) {
                fprintf(stderr, "Error: failed to read intervals for %s in %s\n", chrom, inputs[fi]);
                bmClose(fp);
                goto error_entries;
            }
            for (uint64_t j = 0; j < o->l; j++) {
                if (contextFilterSet && o->context[j] != contextFilter) continue;
                if (o->coverage[j] < minCoverage) continue;
                if (!o->entryid) continue;
                uint32_t id = o->entryid[j];
                if (id == 0) continue;
                if (id >= idmap.n || !idmap.names[id]) continue;
                entries[id].n_sites++;
                entries[id].total_cov += o->coverage[j];
                entries[id].sum_meth_weighted += ((double)o->value[j]) * o->coverage[j];
                if (fp->hdr->version & BM_CONTEXT) {
                    if (o->context[j] == 1) entries[id].n_sites_cg++;
                    else if (o->context[j] == 3) entries[id].n_sites_chh++;
                }
            }
            bmDestroyOverlappingIntervals(o);
        }
        bmClose(fp);
    }

    FILE *out = fopen(outputPath, "w");
    if (!out) {
        fprintf(stderr, "Error: cannot open output file %s: %s\n", outputPath, strerror(errno));
        goto error_entries;
    }
    fprintf(out, "cell_id\tn_sites\ttotal_coverage\tmean_coverage\tmean_meth\tn_sites_cg\tn_sites_chh\tcg_chh_ratio\n");
    for (size_t i = 1; i < idmap.n; i++) {
        if (!idmap.names[i]) continue;
        double mean_cov = entries[i].n_sites ? (double)entries[i].total_cov / (double)entries[i].n_sites : 0.0;
        double mean_meth = entries[i].total_cov ? entries[i].sum_meth_weighted / (double)entries[i].total_cov : 0.0;
        double cg_chh_ratio = 0.0;
        if (entries[i].n_sites_chh > 0) {
            cg_chh_ratio = (double)entries[i].n_sites_cg / (double)entries[i].n_sites_chh;
        } else if (entries[i].n_sites_cg > 0) {
            cg_chh_ratio = (double)entries[i].n_sites_cg;
        }
        fprintf(out, "%s\t%" PRIu64 "\t%" PRIu64 "\t%.6f\t%.6f\t%" PRIu64 "\t%" PRIu64 "\t%.6f\n",
                idmap.names[i],
                entries[i].n_sites,
                entries[i].total_cov,
                mean_cov,
                mean_meth,
                entries[i].n_sites_cg,
                entries[i].n_sites_chh,
                cg_chh_ratio);
    }
    fclose(out);
    goto success;

error_entries:
    dm_sc_qc_free_entries(entries);
    dm_sc_free_idmap(&idmap);
error:
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPath);
cleanup:
    return 1;

success:
    dm_sc_qc_free_entries(entries);
    dm_sc_free_idmap(&idmap);
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPath);
    return 0;
}

typedef struct {
    char *chrom;
    uint32_t start;
    uint32_t end;
    char *name;
    char *id;
} dm_region_t;

typedef struct {
    int row;
    int col;
    uint64_t sum_cov;
    double sum_meth_weighted;
    uint64_t n_sites;
} dm_matrix_accum_t;

typedef struct {
    int row;
    int col;
    double val;
} dm_mtx_entry_t;

typedef struct {
    int row;
    int col;
    double mc;
    double cov;
} dm_count_entry_t;

typedef struct {
    double alpha;
    double beta;
    int valid;
} dm_eb_prior_t;

static double dm_sc_clamp_prob(double p) {
    const double eps = 1e-6;
    if (p < eps) return eps;
    if (p > 1.0 - eps) return 1.0 - eps;
    return p;
}

static double dm_sc_posterior_mean(double alpha, double beta, double k, double n) {
    double denom = alpha + beta + n;
    if (denom == 0.0) return 0.0;
    return (alpha + k) / denom;
}

static double dm_sc_posterior_var(double alpha, double beta, double k, double n) {
    double a = alpha + k;
    double b = beta + (n - k);
    double denom = (alpha + beta + n);
    if (denom == 0.0) return 0.0;
    double denom2 = denom * denom * (denom + 1.0);
    if (denom2 == 0.0) return 0.0;
    return (a * b) / denom2;
}

static void dm_sc_posterior_ci(double alpha, double beta, double k, double n, double *low, double *high) {
    double mean = dm_sc_posterior_mean(alpha, beta, k, n);
    double var = dm_sc_posterior_var(alpha, beta, k, n);
    double sd = var > 0.0 ? sqrt(var) : 0.0;
    double l = mean - 1.96 * sd;
    double h = mean + 1.96 * sd;
    if (l < 0.0) l = 0.0;
    if (h > 1.0) h = 1.0;
    if (low) *low = l;
    if (high) *high = h;
}

typedef struct {
    char *chrom;
    int *idx;
    size_t n;
    size_t cap;
} dm_region_index_t;

typedef struct {
    char *chrom;
    uint32_t binsize;
    uint32_t chrom_len;
    size_t offset;
    size_t nbin;
} dm_bin_index_t;

static void dm_sc_matrix_free_regions(dm_region_t *regions, size_t n) {
    if (!regions) return;
    for (size_t i = 0; i < n; i++) {
        free(regions[i].chrom);
        free(regions[i].name);
        free(regions[i].id);
    }
    free(regions);
}

static void dm_sc_matrix_free_region_index(dm_region_index_t *idx, size_t n) {
    if (!idx) return;
    for (size_t i = 0; i < n; i++) {
        free(idx[i].chrom);
        free(idx[i].idx);
    }
    free(idx);
}

static void dm_sc_matrix_free_bins(dm_bin_index_t *bins, size_t n) {
    if (!bins) return;
    for (size_t i = 0; i < n; i++) free(bins[i].chrom);
    free(bins);
}

static int dm_sc_matrix_find_or_add_id(char ***ids, size_t *n, size_t *cap, const char *id) {
    for (size_t i = 0; i < *n; i++) {
        if (strcmp((*ids)[i], id) == 0) return (int)i;
    }
    if (*n == *cap) {
        size_t newCap = (*cap == 0) ? 16 : (*cap * 2);
        char **tmp = realloc(*ids, newCap * sizeof(char *));
        if (!tmp) return -1;
        *ids = tmp;
        *cap = newCap;
    }
    (*ids)[*n] = strdup(id);
    if (!(*ids)[*n]) return -1;
    (*n)++;
    return (int)(*n - 1);
}

static int dm_sc_matrix_find_or_add_accum(dm_matrix_accum_t **acc, size_t *n, size_t *cap, int row, int col) {
    // Linear search for existing entries; this can be replaced with a hash map for large datasets.
    for (size_t i = 0; i < *n; i++) {
        if ((*acc)[i].row == row && (*acc)[i].col == col) return (int)i;
    }
    if (*n == *cap) {
        size_t newCap = (*cap == 0) ? 64 : (*cap * 2);
        dm_matrix_accum_t *tmp = realloc(*acc, newCap * sizeof(dm_matrix_accum_t));
        if (!tmp) return -1;
        *acc = tmp;
        *cap = newCap;
    }
    (*acc)[*n].row = row;
    (*acc)[*n].col = col;
    (*acc)[*n].sum_cov = 0;
    (*acc)[*n].sum_meth_weighted = 0.0;
    (*acc)[*n].n_sites = 0;
    (*n)++;
    return (int)(*n - 1);
}

static int dm_sc_matrix_region_group(dm_region_index_t **groups, size_t *nGroups, size_t *capGroups, const char *chrom) {
    for (size_t i = 0; i < *nGroups; i++) {
        if (strcmp((*groups)[i].chrom, chrom) == 0) return (int)i;
    }
    if (*nGroups == *capGroups) {
        size_t newCap = (*capGroups == 0) ? 8 : (*capGroups * 2);
        dm_region_index_t *tmp = realloc(*groups, newCap * sizeof(dm_region_index_t));
        if (!tmp) return -1;
        *groups = tmp;
        *capGroups = newCap;
    }
    (*groups)[*nGroups].chrom = strdup(chrom);
    if (!(*groups)[*nGroups].chrom) return -1;
    (*groups)[*nGroups].idx = NULL;
    (*groups)[*nGroups].n = 0;
    (*groups)[*nGroups].cap = 0;
    (*nGroups)++;
    return (int)(*nGroups - 1);
}

static int dm_sc_matrix_group_append(dm_region_index_t *group, int regionIdx) {
    if (group->n == group->cap) {
        size_t newCap = (group->cap == 0) ? 8 : (group->cap * 2);
        int *tmp = realloc(group->idx, newCap * sizeof(int));
        if (!tmp) return -1;
        group->idx = tmp;
        group->cap = newCap;
    }
    group->idx[group->n++] = regionIdx;
    return 0;
}

static int dm_sc_matrix_load_bed(const char *path, dm_region_t **regions, size_t *nRegions, dm_region_index_t **groups, size_t *nGroups) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open BED file %s\n", path);
        return 1;
    }
    size_t capRegions = 0, capGroups = 0;
    char line[4096];
    size_t idxCounter = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *chrom = strtok(line, "\t\n");
        char *startS = strtok(NULL, "\t\n");
        char *endS = strtok(NULL, "\t\n");
        char *nameS = strtok(NULL, "\t\n");
        if (!chrom || !startS || !endS) continue;
        uint32_t start = (uint32_t)strtoul(startS, NULL, 10);
        uint32_t end = (uint32_t)strtoul(endS, NULL, 10);
        if (*nRegions == capRegions) {
            size_t newCap = capRegions ? capRegions * 2 : 64;
            dm_region_t *tmp = realloc(*regions, newCap * sizeof(dm_region_t));
            if (!tmp) { fclose(f); return 1; }
            *regions = tmp; capRegions = newCap;
        }
        dm_region_t *r = &(*regions)[*nRegions];
        r->chrom = strdup(chrom);
        r->start = start;
        r->end = end;
        char idbuf[64];
        snprintf(idbuf, sizeof(idbuf), "region_%06zu", idxCounter + 1);
        r->id = strdup(idbuf);
        r->name = strdup(nameS ? nameS : idbuf);
        if (!r->chrom || !r->id || !r->name) { fclose(f); return 1; }
        int gidx = dm_sc_matrix_region_group(groups, nGroups, &capGroups, chrom);
        if (gidx < 0) { fclose(f); return 1; }
        if (dm_sc_matrix_group_append(&(*groups)[gidx], (int)*nRegions) != 0) { fclose(f); return 1; }
        (*nRegions)++; idxCounter++;
    }
    fclose(f);
    return 0;
}

static dm_bin_index_t *dm_sc_matrix_make_bins(binaMethFile_t *fp, uint32_t binsize, size_t *nBinsOut, dm_region_t **regions, size_t *nRegionsOut) {
    size_t nChrom = fp->cl->nKeys;
    dm_bin_index_t *bins = calloc(nChrom, sizeof(dm_bin_index_t));
    if (!bins) return NULL;
    size_t totalRegions = 0;
    size_t capRegions = 0;
    dm_region_t *reg = NULL;
    char namebuf[128];
    for (size_t i = 0; i < nChrom; i++) {
        bins[i].chrom = strdup(fp->cl->chrom[i]);
        bins[i].binsize = binsize;
        bins[i].chrom_len = fp->cl->len[i];
        bins[i].offset = totalRegions;
        bins[i].nbin = (fp->cl->len[i] + binsize - 1) / binsize;
        for (size_t b = 0; b < bins[i].nbin; b++) {
            if (totalRegions == capRegions) {
                size_t newCap = capRegions ? capRegions * 2 : 128;
                dm_region_t *tmp = realloc(reg, newCap * sizeof(dm_region_t));
                if (!tmp) { dm_sc_matrix_free_bins(bins, nChrom); free(reg); return NULL; }
                reg = tmp; capRegions = newCap;
            }
            dm_region_t *r = &reg[totalRegions];
            r->chrom = strdup(fp->cl->chrom[i]);
            r->start = (uint32_t)(b * binsize);
            uint32_t end = r->start + binsize;
            if (end > fp->cl->len[i]) end = fp->cl->len[i];
            r->end = end;
            snprintf(namebuf, sizeof(namebuf), "%s:%u-%u", fp->cl->chrom[i], r->start, r->end);
            char idbuf[64];
            snprintf(idbuf, sizeof(idbuf), "region_%06zu", totalRegions + 1);
            r->name = strdup(namebuf);
            r->id = strdup(idbuf);
            if (!r->chrom || !r->name || !r->id) { dm_sc_matrix_free_bins(bins, nChrom); dm_sc_matrix_free_regions(reg, totalRegions+1); return NULL; }
            totalRegions++;
        }
    }
    *regions = reg;
    *nRegionsOut = totalRegions;
    *nBinsOut = nChrom;
    return bins;
}

static int dm_sc_matrix_region_lookup_bed(dm_region_index_t *groups, size_t nGroups, dm_region_t *regions, const char *chrom, uint32_t pos) {
    for (size_t i = 0; i < nGroups; i++) {
        if (strcmp(groups[i].chrom, chrom) != 0) continue;
        for (size_t j = 0; j < groups[i].n; j++) {
            int idx = groups[i].idx[j];
            dm_region_t *r = &regions[idx];
            if (pos >= r->start && pos < r->end) return idx;
        }
    }
    return -1;
}

static int dm_sc_matrix_bin_lookup(dm_bin_index_t *bins, size_t nBins, const char *chrom, uint32_t pos, int *colOut) {
    for (size_t i = 0; i < nBins; i++) {
        if (strcmp(bins[i].chrom, chrom) != 0) continue;
        if (pos >= bins[i].chrom_len) return -1;
        size_t binIdx = pos / bins[i].binsize;
        if (binIdx >= bins[i].nbin) return -1;
        *colOut = (int)(bins[i].offset + binIdx);
        return 0;
    }
    return -1;
}

static int dm_sc_matrix_write_features(const char *prefix, dm_region_t *regions, size_t nRegions) {
    size_t len = strlen(prefix) + 16;
    char *path = malloc(len);
    if (!path) return 1;
    snprintf(path, len, "%s.features.tsv", prefix);
    FILE *f = fopen(path, "w");
    free(path);
    if (!f) return 1;
    for (size_t i = 0; i < nRegions; i++) {
        fprintf(f, "%s\t%s\t%u\t%u\t%s\n", regions[i].id, regions[i].chrom, regions[i].start, regions[i].end, regions[i].name);
    }
    fclose(f);
    return 0;
}

static int dm_sc_matrix_write_barcodes(const char *prefix, char **ids, size_t nIds) {
    size_t len = strlen(prefix) + 18;
    char *path = malloc(len);
    if (!path) return 1;
    snprintf(path, len, "%s.barcodes.tsv", prefix);
    FILE *f = fopen(path, "w");
    free(path);
    if (!f) return 1;
    for (size_t i = 0; i < nIds; i++) fprintf(f, "%s\n", ids[i]);
    fclose(f);
    return 0;
}

static int dm_sc_copy_file(const char *src, const char *dst) {
    FILE *in = fopen(src, "r");
    if (!in) return 1;
    FILE *out = fopen(dst, "w");
    if (!out) { fclose(in); return 1; }
    char buf[8192];
    size_t nread = 0;
    while ((nread = fread(buf, 1, sizeof(buf), in)) > 0) {
        if (fwrite(buf, 1, nread, out) != nread) {
            fclose(in);
            fclose(out);
            return 1;
        }
    }
    fclose(in);
    fclose(out);
    return 0;
}

static double dm_sc_matrix_value(dm_matrix_accum_t *a, const char *valueType, const char *agg) {
    if (strcasecmp(valueType, "coverage") == 0) {
        if (strcasecmp(agg, "sum") == 0) return (double)a->sum_cov;
        if (a->n_sites == 0) return 0.0;
        return (double)a->sum_cov / (double)a->n_sites;
    }
    if (a->sum_cov == 0) return 0.0;
    return a->sum_meth_weighted / (double)a->sum_cov;
}

static int dm_mtx_entry_cmp(const void *a, const void *b) {
    const dm_mtx_entry_t *ea = a;
    const dm_mtx_entry_t *eb = b;
    if (ea->row != eb->row) return ea->row < eb->row ? -1 : 1;
    if (ea->col != eb->col) return ea->col < eb->col ? -1 : 1;
    return 0;
}

static int dm_read_mtx(const char *path, int *nRows, int *nCols, dm_mtx_entry_t **entries, size_t *nEntriesOut) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open matrix file %s\n", path);
        return 1;
    }
    char line[4096];
    int gotHeader = 0;
    int rows = 0, cols = 0;
    size_t nnz = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '%') continue;
        if (!gotHeader) {
            if (sscanf(line, "%d %d %zu", &rows, &cols, &nnz) != 3) {
                fprintf(stderr, "Error: invalid matrix header in %s\n", path);
                fclose(f);
                return 1;
            }
            gotHeader = 1;
            break;
        }
    }
    if (!gotHeader) {
        fprintf(stderr, "Error: missing matrix header in %s\n", path);
        fclose(f);
        return 1;
    }
    dm_mtx_entry_t *vals = calloc(nnz, sizeof(dm_mtx_entry_t));
    if (!vals) { fclose(f); return 1; }
    size_t idx = 0;
    while (idx < nnz && fgets(line, sizeof(line), f)) {
        if (line[0] == '%' || line[0] == '\n') continue;
        int r = 0, c = 0;
        double v = 0.0;
        if (sscanf(line, "%d %d %lf", &r, &c, &v) != 3) continue;
        vals[idx].row = r - 1;
        vals[idx].col = c - 1;
        vals[idx].val = v;
        idx++;
    }
    fclose(f);
    if (idx != nnz) {
        fprintf(stderr, "Warning: expected %zu entries in %s, read %zu\n", nnz, path, idx);
        nnz = idx;
    }
    *nRows = rows;
    *nCols = cols;
    *entries = vals;
    *nEntriesOut = nnz;
    return 0;
}


static void dm_sc_matrix_usage() {
    fprintf(stderr,
            "Usage: dmtools sc-matrix -i <input.dm> [-i <input2.dm> ...] -o <output_prefix> (--bed regions.bed | --binsize N) [options]\n"
            "Options:\n"
            "  -i, --input <dm>         Input DM file with IDs (can be repeated)\n"
            "  -o, --output <prefix>    Output prefix for matrix files\n"
            "      --bed <bed>          BED file of regions (chrom start end [name])\n"
            "      --binsize <N>        Fixed bin size for genome-wide bins\n"
            "      --context <CTX>      Context filter: C, CG, CHG, CHH (default: no filter)\n"
            "      --min-coverage <N>   Minimum coverage per site (default: 1)\n"
            "      --value <TYPE>       mean-meth (default) or coverage\n"
            "      --agg <METHOD>       mean (default) or sum (for coverage)\n"
            "      --emit-counts        Emit mC/cov MatrixMarket files for shrinkage\n"
            "      --sparse             Write sparse Matrix Market output (default)\n"
            "      --dense              Write dense TSV matrix\n"
            "      --eb                 Enable empirical Bayes shrinkage (outputs eb_mean + var/CI)\n"
            "      --prior-strength <M> Prior strength M (alpha+beta) for EB (default: 50)\n"
            "  -h, --help               Show this help message\n");
}

int dm_sc_matrix_main(int argc, char **argv) {
    char **inputs = NULL; size_t nInputs = 0, capInputs = 0;
    char *outputPrefix = NULL;
    char *bedPath = NULL;
    uint32_t binsize = 0;
    uint8_t contextFilter = 0; int contextSet = 0;
    uint64_t minCoverage = 1;
    char *valueType = strdup("mean-meth");
    char *aggMethod = strdup("mean");
    int dense = 0; int sparse = 1;
    int ebEnabled = 0;
    int emitCounts = 0;
    double priorStrength = 50.0;

    if (!valueType || !aggMethod) goto error;

    static struct option long_opts[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"bed", required_argument, 0, 1},
        {"binsize", required_argument, 0, 2},
        {"context", required_argument, 0, 3},
        {"min-coverage", required_argument, 0, 4},
        {"value", required_argument, 0, 5},
        {"agg", required_argument, 0, 6},
        {"sparse", no_argument, 0, 7},
        {"dense", no_argument, 0, 8},
        {"eb", no_argument, 0, 9},
        {"emit-counts", no_argument, 0, 10},
        {"prior-strength", required_argument, 0, 11},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int opt, long_idx = 0;
    while ((opt = getopt_long(argc, argv, "hi:o:", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 'i': {
                if (nInputs == capInputs) {
                    size_t newCap = capInputs ? capInputs * 2 : 4;
                    char **tmp = realloc(inputs, newCap * sizeof(char*));
                    if (!tmp) { fprintf(stderr, "Error: memory allocation failed for inputs.\n"); goto error; }
                    inputs = tmp; capInputs = newCap;
                }
                inputs[nInputs] = strdup(optarg);
                if (!inputs[nInputs]) { fprintf(stderr, "Error: memory allocation failed for input path.\n"); goto error; }
                nInputs++;
                break; }
            case 'o':
                outputPrefix = strdup(optarg);
                if (!outputPrefix) { fprintf(stderr, "Error: memory allocation failed for output prefix.\n"); goto error; }
                break;
            case 1:
                bedPath = strdup(optarg);
                if (!bedPath) { fprintf(stderr, "Error: memory allocation failed for bed path.\n"); goto error; }
                break;
            case 2:
                binsize = (uint32_t)strtoul(optarg, NULL, 10);
                break;
            case 3:
                if (dm_sc_qc_parse_context(optarg, &contextFilter) != 0) {
                    fprintf(stderr, "Error: invalid context '%s'. Use C, CG, CHG, or CHH.\n", optarg);
                    goto error;
                }
                contextSet = 1;
                break;
            case 4:
                minCoverage = strtoull(optarg, NULL, 10);
                break;
            case 5:
                free(valueType);
                valueType = strdup(optarg);
                if (!valueType) goto error;
                break;
            case 6:
                free(aggMethod);
                aggMethod = strdup(optarg);
                if (!aggMethod) goto error;
                break;
            case 7:
                sparse = 1; dense = 0;
                break;
            case 8:
                dense = 1; sparse = 0;
                break;
            case 9:
                ebEnabled = 1;
                break;
            case 10:
                emitCounts = 1;
                break;
            case 11:
                priorStrength = atof(optarg);
                if (priorStrength <= 0.0) {
                    fprintf(stderr, "Error: --prior-strength must be positive.\n");
                    goto error;
                }
                break;
            case 'h':
                dm_sc_matrix_usage();
                goto cleanup;
            default:
                dm_sc_matrix_usage();
                goto error;
        }
    }

    if (nInputs == 0 || !outputPrefix) {
        fprintf(stderr, "Error: -i <input.dm> and -o <output_prefix> are required.\n");
        dm_sc_matrix_usage();
        goto error;
    }
    if ((bedPath && binsize > 0) || (!bedPath && binsize == 0)) {
        fprintf(stderr, "Error: specify exactly one of --bed or --binsize.\n");
        dm_sc_matrix_usage();
        goto error;
    }
    if (ebEnabled && strcasecmp(valueType, "mean-meth") != 0) {
        fprintf(stderr, "Warning: --eb only applies to mean-meth; ignoring --value %s.\n", valueType);
    }

    dm_region_t *regions = NULL; size_t nRegions = 0;
    dm_region_index_t *regionGroups = NULL; size_t nGroups = 0;
    dm_bin_index_t *binIndex = NULL; size_t nBinChrom = 0;
    dm_sc_idmap_t idmap = {0};
    char **cellIds = NULL; size_t nCells = 0;
    int *idToRow = NULL;
    uint64_t *cellSites = NULL;
    uint64_t *cellCov = NULL;
    double *cellMeth = NULL;
    dm_matrix_accum_t *acc = NULL; size_t nAcc = 0, capAcc = 0;

    // Prepare regions
    if (bedPath) {
        if (dm_sc_matrix_load_bed(bedPath, &regions, &nRegions, &regionGroups, &nGroups) != 0) {
            fprintf(stderr, "Error: failed to load BED regions.\n");
            goto error_all;
        }
    }

    for (size_t fi = 0; fi < nInputs; fi++) {
        binaMethFile_t *fp = bmOpen(inputs[fi], NULL, "r");
        if (!fp) { fprintf(stderr, "Error: cannot open input file %s\n", inputs[fi]); goto error_all; }
        if (!(fp->hdr->version & BM_ID)) {
            fprintf(stderr, "Error: sc-matrix requires dm files with ID (BM_ID set). File %s has no ID field.\n", inputs[fi]);
            bmClose(fp); goto error_all;
        }
        if (!(fp->hdr->version & BM_COVER)) {
            fprintf(stderr, "Error: sc-matrix requires coverage information. File %s lacks BM_COVER.\n", inputs[fi]);
            bmClose(fp); goto error_all;
        }
        if (contextSet && !(fp->hdr->version & BM_CONTEXT)) {
            fprintf(stderr, "Error: context filter requested but file %s lacks BM_CONTEXT.\n", inputs[fi]);
            bmClose(fp); goto error_all;
        }

        if (!bedPath && !binIndex) {
            binIndex = dm_sc_matrix_make_bins(fp, binsize, &nBinChrom, &regions, &nRegions);
            if (!binIndex) { fprintf(stderr, "Error: failed to generate bins.\n"); bmClose(fp); goto error_all; }
        }

        if (idmap.n == 0) {
            if (dm_sc_load_idmap(inputs[fi], &idmap) != 0) { bmClose(fp); goto error_all; }
            idToRow = calloc(idmap.n, sizeof(int));
            if (!idToRow) { bmClose(fp); goto error_all; }
            for (size_t i = 0; i < idmap.n; i++) idToRow[i] = -1;
            for (size_t i = 1; i < idmap.n; i++) {
                if (!idmap.names[i]) continue;
                idToRow[i] = (int)nCells;
                nCells++;
            }
            cellIds = calloc(nCells, sizeof(char *));
            cellSites = calloc(nCells, sizeof(uint64_t));
            cellCov = calloc(nCells, sizeof(uint64_t));
            cellMeth = calloc(nCells, sizeof(double));
            if (!cellIds || !cellSites || !cellCov || !cellMeth) { bmClose(fp); goto error_all; }
            for (size_t i = 1, row = 0; i < idmap.n; i++) {
                if (!idmap.names[i]) continue;
                cellIds[row++] = strdup(idmap.names[i]);
            }
        }

        for (uint64_t ci = 0; ci < fp->cl->nKeys; ci++) {
            char *chrom = fp->cl->chrom[ci];
            uint32_t chromLen = fp->cl->len[ci];
            bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, chrom, 0, chromLen);
            if (!o) { fprintf(stderr, "Error: failed to read intervals for %s in %s\n", chrom, inputs[fi]); bmClose(fp); goto error_all; }
            for (uint64_t j = 0; j < o->l; j++) {
                if (contextSet && o->context[j] != contextFilter) continue;
                if (o->coverage[j] < minCoverage) continue;
                if (!o->entryid) continue;
                uint32_t id = o->entryid[j];
                if (id == 0) continue;
                if (id >= idmap.n) continue;
                int row = idToRow[id];
                if (row < 0) continue;
                int col = -1;
                if (bedPath) {
                    col = dm_sc_matrix_region_lookup_bed(regionGroups, nGroups, regions, chrom, o->start[j]);
                } else {
                    int tmpCol = -1;
                    if (dm_sc_matrix_bin_lookup(binIndex, nBinChrom, chrom, o->start[j], &tmpCol) == 0) col = tmpCol;
                }
                if (col < 0) continue;
                int accIdx = dm_sc_matrix_find_or_add_accum(&acc, &nAcc, &capAcc, row, col);
                if (accIdx < 0) { bmDestroyOverlappingIntervals(o); bmClose(fp); goto error_all; }
                acc[accIdx].n_sites++;
                acc[accIdx].sum_cov += o->coverage[j];
                acc[accIdx].sum_meth_weighted += ((double)o->value[j]) * o->coverage[j];
                cellSites[row]++;
                cellCov[row] += o->coverage[j];
                cellMeth[row] += ((double)o->value[j]) * o->coverage[j];
            }
            bmDestroyOverlappingIntervals(o);
        }
        bmClose(fp);
    }

    if (!dense && !sparse) sparse = 1;

    double *priorAlpha = NULL;
    double *priorBeta = NULL;
    if (ebEnabled) {
        priorAlpha = calloc(nRegions, sizeof(double));
        priorBeta = calloc(nRegions, sizeof(double));
        if (!priorAlpha || !priorBeta) goto error_all;
        double *sumCov = calloc(nRegions, sizeof(double));
        double *sumMc = calloc(nRegions, sizeof(double));
        if (!sumCov || !sumMc) { free(sumCov); free(sumMc); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            sumCov[acc[i].col] += (double)acc[i].sum_cov;
            sumMc[acc[i].col] += acc[i].sum_meth_weighted;
        }
        for (size_t c = 0; c < nRegions; c++) {
            if (sumCov[c] > 0.0) {
                double pbar = dm_sc_clamp_prob(sumMc[c] / sumCov[c]);
                priorAlpha[c] = pbar * priorStrength;
                priorBeta[c] = (1.0 - pbar) * priorStrength;
            }
        }
        free(sumCov);
        free(sumMc);
    }

    if (dm_sc_matrix_write_features(outputPrefix, regions, nRegions) != 0) {
        fprintf(stderr, "Error: failed to write features.tsv\n");
        goto error_all;
    }
    if (dm_sc_matrix_write_barcodes(outputPrefix, cellIds, nCells) != 0) {
        fprintf(stderr, "Error: failed to write barcodes.tsv\n");
        goto error_all;
    }
    {
        size_t len = strlen(outputPrefix) + 16;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.obs_qc.tsv", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
    fprintf(f, "cell_id\tn_sites\ttotal_coverage\tmean_coverage\tmean_meth\tn_sites_cg\tn_sites_chh\tcg_chh_ratio\n");
    for (size_t i = 0; i < nCells; i++) {
        double mean_cov = cellSites[i] ? (double)cellCov[i] / (double)cellSites[i] : 0.0;
        double mean_meth = cellCov[i] ? cellMeth[i] / (double)cellCov[i] : 0.0;
        fprintf(f, "%s\t%" PRIu64 "\t%" PRIu64 "\t%.6f\t%.6f\t0\t0\t0\n",
                cellIds[i], cellSites[i], cellCov[i], mean_cov, mean_meth);
    }
        fclose(f);
    }

    if (sparse) {
        size_t len = strlen(outputPrefix) + 5;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.matrix.mtx", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            double val = ebEnabled
                ? dm_sc_posterior_mean(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                                       acc[i].sum_meth_weighted, (double)acc[i].sum_cov)
                : dm_sc_matrix_value(&acc[i], valueType, aggMethod);
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, val);
        }
        fclose(f);
    }

    if (emitCounts) {
        size_t len = strlen(outputPrefix) + 9;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.mC.mtx", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, acc[i].sum_meth_weighted);
        }
        fclose(f);

        len = strlen(outputPrefix) + 10;
        path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.cov.mtx", outputPrefix);
        f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, (double)acc[i].sum_cov);
        }
        fclose(f);
    }

    if (sparse && ebEnabled) {
        size_t len = strlen(outputPrefix) + 20;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.matrix.eb_var.mtx", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            double val = dm_sc_posterior_var(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                                             acc[i].sum_meth_weighted, (double)acc[i].sum_cov);
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, val);
        }
        fclose(f);

        len = strlen(outputPrefix) + 22;
        path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.matrix.eb_lower.mtx", outputPrefix);
        f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            double low = 0.0, high = 0.0;
            dm_sc_posterior_ci(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                               acc[i].sum_meth_weighted, (double)acc[i].sum_cov, &low, &high);
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, low);
        }
        fclose(f);

        len = strlen(outputPrefix) + 22;
        path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.matrix.eb_upper.mtx", outputPrefix);
        f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(f, "%zu %zu %zu\n", nCells, nRegions, nAcc);
        for (size_t i = 0; i < nAcc; i++) {
            double low = 0.0, high = 0.0;
            dm_sc_posterior_ci(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                               acc[i].sum_meth_weighted, (double)acc[i].sum_cov, &low, &high);
            fprintf(f, "%d %d %.6f\n", acc[i].row + 1, acc[i].col + 1, high);
        }
        fclose(f);
    }

    if (dense) {
        size_t len = strlen(outputPrefix) + 5;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.tsv", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "cell_id");
        for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%s", regions[c].name);
        fprintf(f, "\n");
        double *matrix = calloc(nCells * nRegions, sizeof(double));
        if (!matrix) { fclose(f); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            double val = ebEnabled
                ? dm_sc_posterior_mean(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                                       acc[i].sum_meth_weighted, (double)acc[i].sum_cov)
                : dm_sc_matrix_value(&acc[i], valueType, aggMethod);
            matrix[((size_t)acc[i].row * nRegions) + acc[i].col] = val;
        }
        for (size_t r = 0; r < nCells; r++) {
            fprintf(f, "%s", cellIds[r]);
            for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%.6f", matrix[r * nRegions + c]);
            fprintf(f, "\n");
        }
        free(matrix);
        fclose(f);
    }

    if (dense && ebEnabled) {
        size_t len = strlen(outputPrefix) + 12;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.eb_var.tsv", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "cell_id");
        for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%s", regions[c].name);
        fprintf(f, "\n");
        double *matrix = calloc(nCells * nRegions, sizeof(double));
        if (!matrix) { fclose(f); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            double val = dm_sc_posterior_var(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                                             acc[i].sum_meth_weighted, (double)acc[i].sum_cov);
            matrix[((size_t)acc[i].row * nRegions) + acc[i].col] = val;
        }
        for (size_t r = 0; r < nCells; r++) {
            fprintf(f, "%s", cellIds[r]);
            for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%.6f", matrix[r * nRegions + c]);
            fprintf(f, "\n");
        }
        free(matrix);
        fclose(f);

        len = strlen(outputPrefix) + 13;
        path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.eb_lower.tsv", outputPrefix);
        f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "cell_id");
        for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%s", regions[c].name);
        fprintf(f, "\n");
        matrix = calloc(nCells * nRegions, sizeof(double));
        if (!matrix) { fclose(f); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            double low = 0.0, high = 0.0;
            dm_sc_posterior_ci(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                               acc[i].sum_meth_weighted, (double)acc[i].sum_cov, &low, &high);
            matrix[((size_t)acc[i].row * nRegions) + acc[i].col] = low;
        }
        for (size_t r = 0; r < nCells; r++) {
            fprintf(f, "%s", cellIds[r]);
            for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%.6f", matrix[r * nRegions + c]);
            fprintf(f, "\n");
        }
        free(matrix);
        fclose(f);

        len = strlen(outputPrefix) + 13;
        path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.eb_upper.tsv", outputPrefix);
        f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "cell_id");
        for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%s", regions[c].name);
        fprintf(f, "\n");
        matrix = calloc(nCells * nRegions, sizeof(double));
        if (!matrix) { fclose(f); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            double low = 0.0, high = 0.0;
            dm_sc_posterior_ci(priorAlpha[acc[i].col], priorBeta[acc[i].col],
                               acc[i].sum_meth_weighted, (double)acc[i].sum_cov, &low, &high);
            matrix[((size_t)acc[i].row * nRegions) + acc[i].col] = high;
        }
        for (size_t r = 0; r < nCells; r++) {
            fprintf(f, "%s", cellIds[r]);
            for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%.6f", matrix[r * nRegions + c]);
            fprintf(f, "\n");
        }
        free(matrix);
        fclose(f);
    }

    goto success;

error_all:
    dm_sc_matrix_free_regions(regions, nRegions);
    dm_sc_matrix_free_region_index(regionGroups, nGroups);
    dm_sc_matrix_free_bins(binIndex, nBinChrom);
    for (size_t i = 0; i < nCells; i++) free(cellIds[i]);
    free(cellIds);
    free(idToRow);
    dm_sc_free_idmap(&idmap);
    free(cellSites);
    free(cellCov);
    free(cellMeth);
    free(acc);
    free(priorAlpha);
    free(priorBeta);
error:
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPrefix); free(bedPath);
    free(valueType); free(aggMethod);
cleanup:
    return 1;

success:
    dm_sc_matrix_free_regions(regions, nRegions);
    dm_sc_matrix_free_region_index(regionGroups, nGroups);
    dm_sc_matrix_free_bins(binIndex, nBinChrom);
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPrefix); free(bedPath);
    for (size_t i = 0; i < nCells; i++) free(cellIds[i]);
    free(cellIds);
    free(idToRow);
    dm_sc_free_idmap(&idmap);
    free(cellSites);
    free(cellCov);
    free(cellMeth);
    free(acc);
    free(priorAlpha);
    free(priorBeta);
    free(valueType); free(aggMethod);
    return 0;
}

typedef struct {
    char *cell;
    int group_idx;
} dm_sc_group_map_t;

typedef struct {
    int group_idx;
    int region_idx;
    uint64_t sum_cov;
    double sum_meth_weighted;
    uint64_t n_sites;
    int *cells;
    size_t n_cells;
    size_t cap_cells;
} dm_sc_group_accum_t;

static int dm_sc_find_or_add_string(char ***ids, size_t *n, size_t *cap, const char *id) {
    for (size_t i = 0; i < *n; i++) {
        if (strcmp((*ids)[i], id) == 0) return (int)i;
    }
    if (*n == *cap) {
        size_t newCap = (*cap == 0) ? 16 : (*cap * 2);
        char **tmp = realloc(*ids, newCap * sizeof(char *));
        if (!tmp) return -1;
        *ids = tmp;
        *cap = newCap;
    }
    (*ids)[*n] = strdup(id);
    if (!(*ids)[*n]) return -1;
    (*n)++;
    return (int)(*n - 1);
}

static int dm_sc_group_add_cell(dm_sc_group_accum_t *acc, int cell_idx) {
    for (size_t i = 0; i < acc->n_cells; i++) {
        if (acc->cells[i] == cell_idx) return 0;
    }
    if (acc->n_cells == acc->cap_cells) {
        size_t newCap = acc->cap_cells ? acc->cap_cells * 2 : 4;
        int *tmp = realloc(acc->cells, newCap * sizeof(int));
        if (!tmp) return -1;
        acc->cells = tmp;
        acc->cap_cells = newCap;
    }
    acc->cells[acc->n_cells++] = cell_idx;
    return 0;
}

static int dm_sc_group_find_or_add_accum(dm_sc_group_accum_t **acc, size_t *n, size_t *cap, int group_idx, int region_idx) {
    for (size_t i = 0; i < *n; i++) {
        if ((*acc)[i].group_idx == group_idx && (*acc)[i].region_idx == region_idx) return (int)i;
    }
    if (*n == *cap) {
        size_t newCap = (*cap == 0) ? 64 : (*cap * 2);
        dm_sc_group_accum_t *tmp = realloc(*acc, newCap * sizeof(dm_sc_group_accum_t));
        if (!tmp) return -1;
        *acc = tmp;
        *cap = newCap;
    }
    (*acc)[*n].group_idx = group_idx;
    (*acc)[*n].region_idx = region_idx;
    (*acc)[*n].sum_cov = 0;
    (*acc)[*n].sum_meth_weighted = 0.0;
    (*acc)[*n].n_sites = 0;
    (*acc)[*n].cells = NULL;
    (*acc)[*n].n_cells = 0;
    (*acc)[*n].cap_cells = 0;
    (*n)++;
    return (int)(*n - 1);
}

static int dm_sc_group_value(dm_sc_group_accum_t *a, const char *valueType, const char *agg, double *out) {
    if (strcasecmp(valueType, "coverage") == 0) {
        if (strcasecmp(agg, "sum") == 0) {
            *out = (double)a->sum_cov;
            return 0;
        }
        if (a->n_sites == 0) {
            *out = 0.0;
            return 0;
        }
        *out = (double)a->sum_cov / (double)a->n_sites;
        return 0;
    }
    if (a->sum_cov == 0) {
        *out = 0.0;
        return 0;
    }
    *out = a->sum_meth_weighted / (double)a->sum_cov;
    return 0;
}

static void dm_sc_group_free_accum(dm_sc_group_accum_t *acc, size_t n) {
    if (!acc) return;
    for (size_t i = 0; i < n; i++) free(acc[i].cells);
    free(acc);
}

static int dm_sc_group_parse_map(const char *path, dm_sc_group_map_t **map, size_t *nMap, char ***groups, size_t *nGroups) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open groups file %s\n", path);
        return 1;
    }
    size_t capMap = 0, capGroups = 0;
    char line[4096];
    int isHeader = 1;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *cell = strtok(line, "\t\n");
        char *group = strtok(NULL, "\t\n");
        if (!cell || !group) continue;
        if (isHeader && (strcasecmp(cell, "cell") == 0 || strcasecmp(cell, "cell_id") == 0 || strcasecmp(cell, "cellid") == 0)) {
            isHeader = 0;
            continue;
        }
        isHeader = 0;
        int gidx = dm_sc_find_or_add_string(groups, nGroups, &capGroups, group);
        if (gidx < 0) { fclose(f); return 1; }
        if (*nMap == capMap) {
            size_t newCap = capMap ? capMap * 2 : 128;
            dm_sc_group_map_t *tmp = realloc(*map, newCap * sizeof(dm_sc_group_map_t));
            if (!tmp) { fclose(f); return 1; }
            *map = tmp; capMap = newCap;
        }
        (*map)[*nMap].cell = strdup(cell);
        (*map)[*nMap].group_idx = gidx;
        if (!(*map)[*nMap].cell) { fclose(f); return 1; }
        (*nMap)++;
    }
    fclose(f);
    if (*nMap == 0) {
        fprintf(stderr, "Error: no mappings read from %s\n", path);
        return 1;
    }
    return 0;
}

static int dm_sc_idmap_find(const dm_sc_idmap_t *map, const char *name) {
    if (!map || !map->names) return -1;
    for (size_t i = 1; i < map->n; i++) {
        if (map->names[i] && strcmp(map->names[i], name) == 0) return (int)i;
    }
    return -1;
}

int dm_sc_aggregate_main(int argc, char **argv) {
    char **inputs = NULL;
    size_t nInputs = 0, capInputs = 0;
    char *outputPrefix = NULL;
    char *bedPath = NULL;
    uint32_t binsize = 0;
    char *groupsPath = NULL;
    int contextSet = 0;
    uint8_t contextFilter = 0;
    uint64_t minCoverage = 1;
    int dense = 0;
    char *valueType = strdup("mean-meth");
    char *aggMethod = strdup("mean");

    static struct option long_opts[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"bed", required_argument, 0, 1},
        {"binsize", required_argument, 0, 2},
        {"groups", required_argument, 0, 3},
        {"context", required_argument, 0, 4},
        {"min-coverage", required_argument, 0, 5},
        {"value", required_argument, 0, 6},
        {"agg", required_argument, 0, 7},
        {"dense", no_argument, 0, 8},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt, long_idx = 0;
    while ((opt = getopt_long(argc, argv, "hi:o:", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 'i': {
                if (nInputs == capInputs) {
                    size_t newCap = capInputs ? capInputs * 2 : 4;
                    char **tmp = realloc(inputs, newCap * sizeof(char *));
                    if (!tmp) { fprintf(stderr, "Error: memory allocation failed for inputs.\n"); goto error; }
                    inputs = tmp; capInputs = newCap;
                }
                inputs[nInputs++] = strdup(optarg);
                if (!inputs[nInputs - 1]) { fprintf(stderr, "Error: memory allocation failed for input path.\n"); goto error; }
                break;
            }
            case 'o':
                outputPrefix = strdup(optarg);
                if (!outputPrefix) { fprintf(stderr, "Error: memory allocation failed for output prefix.\n"); goto error; }
                break;
            case 1:
                bedPath = strdup(optarg);
                if (!bedPath) { fprintf(stderr, "Error: memory allocation failed for bed path.\n"); goto error; }
                break;
            case 2:
                binsize = (uint32_t)strtoul(optarg, NULL, 10);
                break;
            case 3:
                groupsPath = strdup(optarg);
                if (!groupsPath) { fprintf(stderr, "Error: memory allocation failed for groups path.\n"); goto error; }
                break;
            case 4:
                if (dm_sc_qc_parse_context(optarg, &contextFilter) != 0) {
                    fprintf(stderr, "Error: invalid context '%s'. Use C, CG, CHG, or CHH.\n", optarg);
                    goto error;
                }
                contextSet = 1;
                break;
            case 5:
                minCoverage = strtoull(optarg, NULL, 10);
                break;
            case 6:
                free(valueType);
                valueType = strdup(optarg);
                if (!valueType) goto error;
                break;
            case 7:
                free(aggMethod);
                aggMethod = strdup(optarg);
                if (!aggMethod) goto error;
                break;
            case 8:
                dense = 1;
                break;
            case 'h':
                fprintf(stderr, "%s\n", Help_String_scaggregate);
                goto success;
            default:
                fprintf(stderr, "%s\n", Help_String_scaggregate);
                goto error;
        }
    }

    if (nInputs == 0 || !outputPrefix || !groupsPath || ((bedPath == NULL) == (binsize == 0))) {
        fprintf(stderr, "Error: -i, -o, --groups and exactly one of --bed or --binsize are required.\n");
        fprintf(stderr, "%s\n", Help_String_scaggregate);
        goto error;
    }

    if (strcasecmp(valueType, "mean-meth") != 0 && strcasecmp(valueType, "coverage") != 0) {
        fprintf(stderr, "Error: --value must be mean-meth or coverage.\n");
        goto error;
    }
    if (strcasecmp(aggMethod, "mean") != 0 && strcasecmp(aggMethod, "sum") != 0) {
        fprintf(stderr, "Error: --agg must be mean or sum.\n");
        goto error;
    }

    dm_sc_group_map_t *mapping = NULL;
    size_t nMap = 0;
    char **groups = NULL;
    size_t nGroups = 0, capGroups = 0;
    if (dm_sc_group_parse_map(groupsPath, &mapping, &nMap, &groups, &nGroups) != 0) goto error;

    dm_region_t *regions = NULL;
    size_t nRegions = 0;
    dm_region_index_t *regionGroups = NULL;
    size_t nRegionGroups = 0;
    dm_bin_index_t *binIndex = NULL;
    size_t nBinChrom = 0;

    dm_sc_group_accum_t *acc = NULL;
    size_t nAcc = 0, capAcc = 0;
    dm_sc_idmap_t idmap = {0};
    char **cellIds = NULL; size_t nCells = 0;
    int *idToRow = NULL;
    int *idToGroup = NULL;

    for (size_t fi = 0; fi < nInputs; fi++) {
        binaMethFile_t *fp = bmOpen(inputs[fi], NULL, "r");
        if (!fp) { fprintf(stderr, "Error: cannot open input file %s\n", inputs[fi]); goto error_all; }
        if (!(fp->hdr->version & BM_ID)) { fprintf(stderr, "Error: sc-aggregate requires dm files with ID (BM_ID set). File %s has no ID field.\n", inputs[fi]); bmClose(fp); goto error_all; }
        if (!(fp->hdr->version & BM_COVER)) { fprintf(stderr, "Error: sc-aggregate requires coverage information. File %s lacks BM_COVER.\n", inputs[fi]); bmClose(fp); goto error_all; }
        if (contextSet && !(fp->hdr->version & BM_CONTEXT)) { fprintf(stderr, "Error: context filter requested but file %s lacks BM_CONTEXT.\n", inputs[fi]); bmClose(fp); goto error_all; }

        if (idmap.n == 0) {
            if (dm_sc_load_idmap(inputs[fi], &idmap) != 0) { bmClose(fp); goto error_all; }
            idToRow = calloc(idmap.n, sizeof(int));
            idToGroup = calloc(idmap.n, sizeof(int));
            if (!idToRow || !idToGroup) { bmClose(fp); goto error_all; }
            for (size_t i = 0; i < idmap.n; i++) { idToRow[i] = -1; idToGroup[i] = -1; }
            for (size_t i = 1; i < idmap.n; i++) {
                if (!idmap.names[i]) continue;
                idToRow[i] = (int)nCells;
                nCells++;
            }
            cellIds = calloc(nCells, sizeof(char *));
            if (!cellIds) { bmClose(fp); goto error_all; }
            for (size_t i = 1, row = 0; i < idmap.n; i++) {
                if (!idmap.names[i]) continue;
                cellIds[row++] = strdup(idmap.names[i]);
            }
            for (size_t i = 0; i < nMap; i++) {
                int id = dm_sc_idmap_find(&idmap, mapping[i].cell);
                if (id > 0 && (size_t)id < idmap.n) {
                    idToGroup[id] = mapping[i].group_idx;
                }
            }
        }

        if (fi == 0) {
            if (bedPath) {
                if (dm_sc_matrix_load_bed(bedPath, &regions, &nRegions, &regionGroups, &nRegionGroups) != 0) { bmClose(fp); goto error_all; }
            } else {
                binIndex = dm_sc_matrix_make_bins(fp, binsize, &nBinChrom, &regions, &nRegions);
                if (!binIndex) { fprintf(stderr, "Error: failed to generate bins.\n"); bmClose(fp); goto error_all; }
            }
        }

        int missing_warned = 0;
        for (uint64_t ci = 0; ci < fp->cl->nKeys; ci++) {
            char *chrom = fp->cl->chrom[ci];
            uint32_t chromLen = fp->cl->len[ci];
            bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, chrom, 0, chromLen);
            if (!o) { fprintf(stderr, "Error: failed to read intervals for %s in %s\n", chrom, inputs[fi]); bmClose(fp); goto error_all; }
            for (uint64_t j = 0; j < o->l; j++) {
                if (contextSet && o->context[j] != contextFilter) continue;
                if (o->coverage[j] < minCoverage) continue;
                if (!o->entryid) continue;
                uint32_t id = o->entryid[j];
                if (id == 0) continue;
                if (id >= idmap.n) continue;
                int gidx = idToGroup[id];
                if (gidx < 0) {
                    if (!missing_warned) {
                        fprintf(stderr, "Warning: cell_id %u not found in groups mapping, skipping (further warnings suppressed)\n", id);
                        missing_warned = 1;
                    }
                    continue;
                }
                int rowCell = idToRow[id];
                if (rowCell < 0) { bmDestroyOverlappingIntervals(o); bmClose(fp); goto error_all; }
                int col = -1;
                if (bedPath) {
                    col = dm_sc_matrix_region_lookup_bed(regionGroups, nRegionGroups, regions, chrom, o->start[j]);
                } else {
                    int tmpCol = -1;
                    if (dm_sc_matrix_bin_lookup(binIndex, nBinChrom, chrom, o->start[j], &tmpCol) == 0) col = tmpCol;
                }
                if (col < 0) continue;
                int accIdx = dm_sc_group_find_or_add_accum(&acc, &nAcc, &capAcc, gidx, col);
                if (accIdx < 0) { bmDestroyOverlappingIntervals(o); bmClose(fp); goto error_all; }
                acc[accIdx].n_sites++;
                acc[accIdx].sum_cov += o->coverage[j];
                acc[accIdx].sum_meth_weighted += ((double)o->value[j]) * o->coverage[j];
                if (dm_sc_group_add_cell(&acc[accIdx], rowCell) != 0) { bmDestroyOverlappingIntervals(o); bmClose(fp); goto error_all; }
            }
            bmDestroyOverlappingIntervals(o);
        }
        bmClose(fp);
    }

    size_t pathLen = strlen(outputPrefix) + 13;
    char *groupPath = malloc(pathLen);
    if (!groupPath) goto error_all;
    snprintf(groupPath, pathLen, "%s.groups.tsv", outputPrefix);
    FILE *gout = fopen(groupPath, "w");
    free(groupPath);
    if (!gout) { fprintf(stderr, "Error: cannot open output file %s.groups.tsv\n", outputPrefix); goto error_all; }
    fprintf(gout, "group\tchrom\tstart\tend\tn_cells\tn_sites\ttotal_coverage\tmean_coverage\tmean_meth\n");
    for (size_t i = 0; i < nAcc; i++) {
        dm_sc_group_accum_t *a = &acc[i];
        dm_region_t *r = &regions[a->region_idx];
        double value = 0.0;
        dm_sc_group_value(a, valueType, aggMethod, &value);
        double mean_cov = a->n_sites ? (double)a->sum_cov / (double)a->n_sites : 0.0;
        fprintf(gout, "%s\t%s\t%u\t%u\t%zu\t%" PRIu64 "\t%" PRIu64 "\t%.6f\t%.6f\n",
                groups[a->group_idx], r->chrom, r->start, r->end,
                a->n_cells, a->n_sites, a->sum_cov, mean_cov, value);
    }
    fclose(gout);

    if (dense) {
        size_t len = strlen(outputPrefix) + 14;
        char *path = malloc(len);
        if (!path) goto error_all;
        snprintf(path, len, "%s.matrix.tsv", outputPrefix);
        FILE *f = fopen(path, "w");
        free(path);
        if (!f) goto error_all;
        fprintf(f, "group");
        for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%s", regions[c].name);
        fprintf(f, "\n");
        double *matrix = calloc(nGroups * nRegions, sizeof(double));
        if (!matrix) { fclose(f); goto error_all; }
        for (size_t i = 0; i < nAcc; i++) {
            double val = 0.0;
            dm_sc_group_value(&acc[i], valueType, aggMethod, &val);
            matrix[((size_t)acc[i].group_idx * nRegions) + acc[i].region_idx] = val;
        }
        for (size_t g = 0; g < nGroups; g++) {
            fprintf(f, "%s", groups[g]);
            for (size_t c = 0; c < nRegions; c++) fprintf(f, "\t%.6f", matrix[g * nRegions + c]);
            fprintf(f, "\n");
        }
        free(matrix);
        fclose(f);
        if (dm_sc_matrix_write_features(outputPrefix, regions, nRegions) != 0) { fprintf(stderr, "Error: failed to write features.tsv\n"); goto error_all; }
        size_t groupsOnlyLen = strlen(outputPrefix) + 18;
        char *groupsOnly = malloc(groupsOnlyLen);
        if (!groupsOnly) goto error_all;
        snprintf(groupsOnly, groupsOnlyLen, "%s.groups_only.tsv", outputPrefix);
        FILE *gf = fopen(groupsOnly, "w");
        free(groupsOnly);
        if (!gf) goto error_all;
        for (size_t g = 0; g < nGroups; g++) fprintf(gf, "%s\n", groups[g]);
        fclose(gf);
    }

    goto success;

error_all:
    dm_sc_matrix_free_regions(regions, nRegions);
    dm_sc_matrix_free_region_index(regionGroups, nRegionGroups);
    dm_sc_matrix_free_bins(binIndex, nBinChrom);
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPrefix); free(bedPath); free(groupsPath);
    free(valueType); free(aggMethod);
    for (size_t i = 0; i < nMap; i++) free(mapping[i].cell);
    free(mapping);
    for (size_t i = 0; i < nGroups; i++) free(groups[i]);
    free(groups);
    for (size_t i = 0; i < nCells; i++) free(cellIds[i]);
    free(cellIds);
    free(idToRow);
    free(idToGroup);
    dm_sc_free_idmap(&idmap);
    dm_sc_group_free_accum(acc, nAcc);
    goto cleanup;

error:
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPrefix); free(bedPath); free(groupsPath);
    free(valueType); free(aggMethod);
cleanup:
    return 1;

success:
    dm_sc_matrix_free_regions(regions, nRegions);
    dm_sc_matrix_free_region_index(regionGroups, nRegionGroups);
    dm_sc_matrix_free_bins(binIndex, nBinChrom);
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPrefix); free(bedPath); free(groupsPath);
    free(valueType); free(aggMethod);
    for (size_t i = 0; i < nMap; i++) free(mapping[i].cell);
    free(mapping);
    for (size_t i = 0; i < nGroups; i++) free(groups[i]);
    free(groups);
    for (size_t i = 0; i < nCells; i++) free(cellIds[i]);
    free(cellIds);
    free(idToRow);
    free(idToGroup);
    dm_sc_free_idmap(&idmap);
    dm_sc_group_free_accum(acc, nAcc);
    return 0;
}

static void dm_sc_export_usage() {
    fprintf(stderr,
            "Usage: dmtools sc-export -i <input.dm> -o <output_prefix> (--bed regions.bed | --binsize N) [--context CG] [--min-coverage 1] [--to h5ad] [--h5ad-script path]\n"
            "Notes:\n"
            "  - Outputs MatrixMarket bundle (matrix.mtx, barcodes.tsv, features.tsv, obs_qc.tsv).\n"
            "  - --to h5ad calls scripts/sc_export_h5ad.py (requires python3, anndata, scipy, pandas).\n");
}

static int dm_sc_resolve_h5ad_script(const char *argv0, const char *override, char *outPath, size_t outSize) {
    if (!outPath || outSize == 0) return 1;
    outPath[0] = '\0';
    if (override && override[0] != '\0') {
        snprintf(outPath, outSize, "%s", override);
        return access(outPath, R_OK) == 0 ? 0 : 1;
    }

    char exePath[PATH_MAX];
    ssize_t n = readlink("/proc/self/exe", exePath, sizeof(exePath) - 1);
    if (n > 0) {
        exePath[n] = '\0';
        char *slash = strrchr(exePath, '/');
        if (slash) {
            *slash = '\0';
            snprintf(outPath, outSize, "%s/scripts/sc_export_h5ad.py", exePath);
            if (access(outPath, R_OK) == 0) return 0;
        }
    }

    if (argv0 && argv0[0] != '\0' && strchr(argv0, '/')) {
        char resolved[PATH_MAX];
        if (realpath(argv0, resolved)) {
            char *slash = strrchr(resolved, '/');
            if (slash) {
                *slash = '\0';
                snprintf(outPath, outSize, "%s/scripts/sc_export_h5ad.py", resolved);
                if (access(outPath, R_OK) == 0) return 0;
            }
        }
    }

    snprintf(outPath, outSize, "scripts/sc_export_h5ad.py");
    return access(outPath, R_OK) == 0 ? 0 : 1;
}

int dm_sc_export_main(int argc, char **argv) {
    char *input = NULL;
    char *outputPrefix = NULL;
    char *bedPath = NULL;
    char *binsize = NULL;
    char *context = NULL;
    char *minCov = NULL;
    char *toFmt = NULL;
    char *h5adScript = NULL;

    static struct option long_opts[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"bed", required_argument, 0, 1},
        {"binsize", required_argument, 0, 2},
        {"context", required_argument, 0, 3},
        {"min-coverage", required_argument, 0, 4},
        {"to", required_argument, 0, 5},
        {"h5ad-script", required_argument, 0, 6},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt, long_idx = 0;
    int rc = 0;
    while ((opt = getopt_long(argc, argv, "hi:o:", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 'i': input = strdup(optarg); break;
            case 'o': outputPrefix = strdup(optarg); break;
            case 1: bedPath = strdup(optarg); break;
            case 2: binsize = strdup(optarg); break;
            case 3: context = strdup(optarg); break;
            case 4: minCov = strdup(optarg); break;
            case 5: toFmt = strdup(optarg); break;
            case 6: h5adScript = strdup(optarg); break;
            case 'h': dm_sc_export_usage(); rc = 0; goto cleanup;
            default: dm_sc_export_usage(); rc = 1; goto cleanup;
        }
    }

    if (!input || !outputPrefix || ((bedPath == NULL) == (binsize == NULL))) {
        dm_sc_export_usage();
        rc = 1;
        goto cleanup;
    }

    char *argv_sc[32];
    int idx = 0;
    argv_sc[idx++] = "sc-matrix";
    argv_sc[idx++] = "-i"; argv_sc[idx++] = input;
    argv_sc[idx++] = "-o"; argv_sc[idx++] = outputPrefix;
    if (bedPath) { argv_sc[idx++] = "--bed"; argv_sc[idx++] = bedPath; }
    if (binsize) { argv_sc[idx++] = "--binsize"; argv_sc[idx++] = binsize; }
    if (context) { argv_sc[idx++] = "--context"; argv_sc[idx++] = context; }
    if (minCov) { argv_sc[idx++] = "--min-coverage"; argv_sc[idx++] = minCov; }
    argv_sc[idx++] = "--sparse";
    argv_sc[idx++] = "--emit-counts";
    argv_sc[idx++] = "--eb";
    argv_sc[idx] = NULL;

    optind = 1;
    opterr = 1;
    optopt = 0;
    rc = dm_sc_matrix_main(idx, argv_sc);
    if (rc != 0) goto cleanup;

    if (toFmt && strcasecmp(toFmt, "h5ad") == 0) {
        char scriptPath[PATH_MAX];
        if (dm_sc_resolve_h5ad_script(argv[0], h5adScript, scriptPath, sizeof(scriptPath)) != 0) {
            fprintf(stderr, "Error: unable to locate sc_export_h5ad.py. Provide --h5ad-script.\n");
            rc = 1;
            goto cleanup;
        }
        char cmd[PATH_MAX * 2];
        snprintf(cmd, sizeof(cmd),
                 "python3 %s --prefix %s --output %s.h5ad",
                 scriptPath, outputPrefix, outputPrefix);
        int sysrc = system(cmd);
        if (sysrc != 0) {
            int exitCode = -1;
            if (sysrc != -1 && WIFEXITED(sysrc)) exitCode = WEXITSTATUS(sysrc);
            if (exitCode == 2 || exitCode == 127) {
                fprintf(stderr, "Warning: sc-export --to h5ad skipped. Install python3, anndata, scipy, pandas to enable h5ad export.\n");
            } else {
                fprintf(stderr, "Error: sc-export --to h5ad failed (exit code %d).\n", exitCode);
                rc = 1;
                goto cleanup;
            }
        }
    }
cleanup:
    free(input);
    free(outputPrefix);
    free(bedPath);
    free(binsize);
    free(context);
    free(minCov);
    free(toFmt);
    free(h5adScript);
    return rc;
}

typedef struct {
    uint32_t *start;
    uint32_t *end;
    uint16_t *cov;
    double *meth;
    uint8_t *strand;
    uint8_t *context;
    size_t n;
    size_t cap;
    uint32_t last_pos;
} dm_pb_buf_t;

static int dm_pb_buf_append(dm_pb_buf_t *buf, uint32_t pos, uint16_t cov, double meth, uint8_t strand, uint8_t context) {
    if (buf->n > 0 && buf->last_pos == pos) {
        buf->cov[buf->n - 1] += cov;
        buf->meth[buf->n - 1] += meth;
        return 0;
    }
    if (buf->n == buf->cap) {
        size_t newCap = buf->cap ? buf->cap * 2 : 1024;
        uint32_t *ns = realloc(buf->start, newCap * sizeof(uint32_t));
        uint32_t *ne = realloc(buf->end, newCap * sizeof(uint32_t));
        uint16_t *nc = realloc(buf->cov, newCap * sizeof(uint16_t));
        double *nm = realloc(buf->meth, newCap * sizeof(double));
        uint8_t *nstrand = realloc(buf->strand, newCap * sizeof(uint8_t));
        uint8_t *ncontext = realloc(buf->context, newCap * sizeof(uint8_t));
        if (!ns || !ne || !nc || !nm || !nstrand || !ncontext) return 1;
        buf->start = ns; buf->end = ne; buf->cov = nc; buf->meth = nm; buf->strand = nstrand; buf->context = ncontext;
        buf->cap = newCap;
    }
    buf->start[buf->n] = pos;
    buf->end[buf->n] = pos + 1;
    buf->cov[buf->n] = cov;
    buf->meth[buf->n] = meth;
    buf->strand[buf->n] = strand;
    buf->context[buf->n] = context;
    buf->last_pos = pos;
    buf->n++;
    return 0;
}

static void dm_pb_buf_reset(dm_pb_buf_t *buf) {
    buf->n = 0;
    buf->last_pos = (uint32_t)-1;
}

static void dm_pb_buf_free(dm_pb_buf_t *buf) {
    free(buf->start);
    free(buf->end);
    free(buf->cov);
    free(buf->meth);
    free(buf->strand);
    free(buf->context);
}

static void dm_sc_pseudobulk_usage() {
    fprintf(stderr,
            "Usage: dmtools sc-pseudobulk -i <input.dm> -o <out_prefix> --groups <mapping.tsv> [--context CG] [--format dm|bigwig|both] [--coverage]\n");
}

int dm_sc_pseudobulk_main(int argc, char **argv) {
    char *input = NULL;
    char *outputPrefix = NULL;
    char *groupsPath = NULL;
    int contextSet = 0;
    uint8_t contextFilter = 0;
    char *format = strdup("both");
    int writeCoverageTrack = 1;

    static struct option long_opts[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"groups", required_argument, 0, 1},
        {"context", required_argument, 0, 2},
        {"format", required_argument, 0, 3},
        {"coverage", no_argument, 0, 4},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt, long_idx = 0;
    while ((opt = getopt_long(argc, argv, "hi:o:", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 'i': input = strdup(optarg); break;
            case 'o': outputPrefix = strdup(optarg); break;
            case 1: groupsPath = strdup(optarg); break;
            case 2:
                if (dm_sc_qc_parse_context(optarg, &contextFilter) != 0) {
                    fprintf(stderr, "Error: invalid context '%s'. Use C, CG, CHG, or CHH.\n", optarg);
                    return 1;
                }
                contextSet = 1;
                break;
            case 3:
                free(format);
                format = strdup(optarg);
                break;
            case 4:
                writeCoverageTrack = 1;
                break;
            case 'h': dm_sc_pseudobulk_usage(); return 0;
            default: dm_sc_pseudobulk_usage(); return 1;
        }
    }

    if (!input || !outputPrefix || !groupsPath) {
        dm_sc_pseudobulk_usage();
        return 1;
    }

    dm_sc_group_map_t *mapping = NULL;
    size_t nMap = 0;
    char **groups = NULL;
    size_t nGroups = 0, capGroups = 0;
    if (dm_sc_group_parse_map(groupsPath, &mapping, &nMap, &groups, &nGroups) != 0) return 1;

    binaMethFile_t *fp = bmOpen(input, NULL, "r");
    if (!fp) { fprintf(stderr, "Error: cannot open input %s\n", input); return 1; }
    if (!(fp->hdr->version & BM_ID)) {
        fprintf(stderr, "Error: sc-pseudobulk requires DM with ID field.\n");
        bmClose(fp);
        return 1;
    }

    dm_sc_idmap_t idmap = {0};
    if (dm_sc_load_idmap(input, &idmap) != 0) { bmClose(fp); return 1; }
    int *idToGroup = calloc(idmap.n, sizeof(int));
    if (!idToGroup) { bmClose(fp); dm_sc_free_idmap(&idmap); return 1; }
    for (size_t i = 0; i < idmap.n; i++) idToGroup[i] = -1;
    for (size_t i = 0; i < nMap; i++) {
        int id = dm_sc_idmap_find(&idmap, mapping[i].cell);
        if (id > 0 && (size_t)id < idmap.n) idToGroup[id] = mapping[i].group_idx;
    }

    binaMethFile_t **outs = calloc(nGroups, sizeof(binaMethFile_t*));
    dm_pb_buf_t *buffers = calloc(nGroups, sizeof(dm_pb_buf_t));
    if (!outs || !buffers) { bmClose(fp); dm_sc_free_idmap(&idmap); free(idToGroup); free(format); return 1; }

    uint32_t write_type = BM_COVER | BM_END;
    if (fp->hdr->version & BM_STRAND) write_type |= BM_STRAND;
    if (fp->hdr->version & BM_CONTEXT) write_type |= BM_CONTEXT;

    int emitDM = (strcasecmp(format, "dm") == 0) || (strcasecmp(format, "both") == 0);
    int emitBW = (strcasecmp(format, "bigwig") == 0) || (strcasecmp(format, "both") == 0);

    if (emitDM) {
        for (size_t g = 0; g < nGroups; g++) {
            char path[PATH_MAX];
            snprintf(path, sizeof(path), "%s.%s.dm", outputPrefix, groups[g]);
            outs[g] = bmOpen(path, NULL, "w");
            if (!outs[g]) { fprintf(stderr, "Error: cannot open output %s\n", path); goto error; }
            outs[g]->type = write_type;
            if (bmCreateHdr(outs[g], 0)) { fprintf(stderr, "Error: bmCreateHdr failed\n"); goto error; }
            outs[g]->cl = bmCreateChromList(fp->cl->chrom, fp->cl->len, fp->cl->nKeys);
            if (!outs[g]->cl) { fprintf(stderr, "Error: bmCreateChromList failed\n"); goto error; }
            if (bmWriteHdr(outs[g])) { fprintf(stderr, "Error: bmWriteHdr failed\n"); goto error; }
        }
    }

    binaMethFile_t **bw_ml = NULL;
    binaMethFile_t **bw_cov = NULL;
    if (emitBW) {
        bw_ml = calloc(nGroups, sizeof(binaMethFile_t*));
        bw_cov = calloc(nGroups, sizeof(binaMethFile_t*));
        if (!bw_ml || !bw_cov) { goto error; }
        for (size_t g = 0; g < nGroups; g++) {
            char mlPath[PATH_MAX];
            snprintf(mlPath, sizeof(mlPath), "%s.%s.ml.bw", outputPrefix, groups[g]);
            bw_ml[g] = bmOpen(mlPath, NULL, "w");
            if (!bw_ml[g]) { fprintf(stderr, "Error: cannot open %s\n", mlPath); goto error; }
            bw_ml[g]->type = BM_END;
            if (bmCreateHdr(bw_ml[g], 0)) { fprintf(stderr, "Error: bmCreateHdr failed for %s\n", mlPath); goto error; }
            bw_ml[g]->cl = bmCreateChromList(fp->cl->chrom, fp->cl->len, fp->cl->nKeys);
            if (!bw_ml[g]->cl) { fprintf(stderr, "Error: bmCreateChromList failed for %s\n", mlPath); goto error; }
            if (bmWriteHdr(bw_ml[g])) { fprintf(stderr, "Error: bmWriteHdr failed for %s\n", mlPath); goto error; }

            if (writeCoverageTrack && (fp->hdr->version & BM_COVER)) {
                char covPath[PATH_MAX];
                snprintf(covPath, sizeof(covPath), "%s.%s.cov.bw", outputPrefix, groups[g]);
                bw_cov[g] = bmOpen(covPath, NULL, "w");
                if (!bw_cov[g]) { fprintf(stderr, "Error: cannot open %s\n", covPath); goto error; }
                bw_cov[g]->type = BM_END;
                if (bmCreateHdr(bw_cov[g], 0)) { fprintf(stderr, "Error: bmCreateHdr failed for %s\n", covPath); goto error; }
                bw_cov[g]->cl = bmCreateChromList(fp->cl->chrom, fp->cl->len, fp->cl->nKeys);
                if (!bw_cov[g]->cl) { fprintf(stderr, "Error: bmCreateChromList failed for %s\n", covPath); goto error; }
                if (bmWriteHdr(bw_cov[g])) { fprintf(stderr, "Error: bmWriteHdr failed for %s\n", covPath); goto error; }
            }
        }
    }

    for (uint64_t ci = 0; ci < fp->cl->nKeys; ci++) {
        char *chrom = fp->cl->chrom[ci];
        uint32_t chromLen = fp->cl->len[ci];
        bmOverlappingIntervals_t *o = bmGetOverlappingIntervals(fp, chrom, 0, chromLen);
        if (!o) { fprintf(stderr, "Error: failed to read intervals for %s\n", chrom); goto error; }
        for (size_t g = 0; g < nGroups; g++) dm_pb_buf_reset(&buffers[g]);
        for (uint64_t j = 0; j < o->l; j++) {
            if (contextSet && o->context[j] != contextFilter) continue;
            if (!o->entryid) continue;
            uint32_t id = o->entryid[j];
            if (id == 0) continue;
            if (id >= idmap.n) continue;
            int gidx = idToGroup[id];
            if (gidx < 0) continue;
            uint16_t cov = (fp->hdr->version & BM_COVER) ? o->coverage[j] : 1;
            double meth = ((double)o->value[j]) * cov;
            uint8_t strand = (fp->hdr->version & BM_STRAND) ? o->strand[j] : 0;
            uint8_t context = (fp->hdr->version & BM_CONTEXT) ? o->context[j] : 0;
            if (dm_pb_buf_append(&buffers[gidx], o->start[j], cov, meth, strand, context) != 0) {
                bmDestroyOverlappingIntervals(o);
                goto error;
            }
        }
        bmDestroyOverlappingIntervals(o);
        for (size_t g = 0; g < nGroups; g++) {
            if (buffers[g].n == 0) continue;
            float *vals = calloc(buffers[g].n, sizeof(float));
            if (!vals) goto error;
            for (size_t i = 0; i < buffers[g].n; i++) {
                if (buffers[g].cov[i] == 0) vals[i] = 0.0f;
                else vals[i] = (float)(buffers[g].meth[i] / (double)buffers[g].cov[i]);
            }
            char **chroms = calloc(buffers[g].n, sizeof(char*));
            if (!chroms) { free(vals); goto error; }
            for (size_t i = 0; i < buffers[g].n; i++) chroms[i] = chrom;
            if (emitDM) {
                int rc = bmAddIntervals(outs[g], chroms, buffers[g].start, buffers[g].end, vals, buffers[g].cov,
                                        buffers[g].strand, buffers[g].context, NULL, (uint32_t)buffers[g].n);
                if (rc != 0) { fprintf(stderr, "Error: bmAddIntervals failed for group %s\n", groups[g]); free(vals); free(chroms); goto error; }
            }
            if (emitBW && bw_ml[g]) {
                uint32_t *starts_bw = calloc(buffers[g].n, sizeof(uint32_t));
                uint32_t *ends_bw = calloc(buffers[g].n, sizeof(uint32_t));
                if (!starts_bw || !ends_bw) { free(vals); free(chroms); free(starts_bw); free(ends_bw); goto error; }
                for (size_t i = 0; i < buffers[g].n; i++) {
                    starts_bw[i] = (buffers[g].start[i] > 0) ? buffers[g].start[i] - 1 : 0;
                    ends_bw[i] = (buffers[g].end[i] > 0) ? buffers[g].end[i] - 1 : 0;
                }
                if (bmAddIntervals(bw_ml[g], chroms, starts_bw, ends_bw, vals, NULL, NULL, NULL, NULL, (uint32_t)buffers[g].n) != 0) {
                    fprintf(stderr, "Error: failed to write bigWig ml for group %s\n", groups[g]); free(starts_bw); free(ends_bw); free(vals); free(chroms); goto error;
                }
                if (bw_cov[g]) {
                    float *covVals = calloc(buffers[g].n, sizeof(float));
                    if (!covVals) { free(starts_bw); free(ends_bw); free(vals); free(chroms); goto error; }
                    for (size_t i = 0; i < buffers[g].n; i++) covVals[i] = (float)buffers[g].cov[i];
                    if (bmAddIntervals(bw_cov[g], chroms, starts_bw, ends_bw, covVals, NULL, NULL, NULL, NULL, (uint32_t)buffers[g].n) != 0) {
                        fprintf(stderr, "Error: failed to write bigWig cov for group %s\n", groups[g]); free(covVals); free(starts_bw); free(ends_bw); free(vals); free(chroms); goto error;
                    }
                    free(covVals);
                }
                free(starts_bw);
                free(ends_bw);
            }
            free(vals);
            free(chroms);
        }
    }

    bmClose(fp);
    for (size_t g = 0; g < nGroups; g++) {
        if (outs && outs[g]) bmClose(outs[g]);
        if (bw_ml && bw_ml[g]) bmClose(bw_ml[g]);
        if (bw_cov && bw_cov[g]) bmClose(bw_cov[g]);
        dm_pb_buf_free(&buffers[g]);
    }
    dm_sc_free_idmap(&idmap);
    free(idToGroup);
    free(outs);
    free(buffers);
    free(bw_ml);
    free(bw_cov);
    for (size_t i = 0; i < nMap; i++) free(mapping[i].cell);
    free(mapping);
    for (size_t i = 0; i < nGroups; i++) free(groups[i]);
    free(groups);
    free(input);
    free(outputPrefix);
    free(groupsPath);
    free(format);
    return 0;

error:
    bmClose(fp);
    for (size_t g = 0; g < nGroups; g++) {
        if (outs && outs[g]) bmClose(outs[g]);
        if (bw_ml && bw_ml[g]) bmClose(bw_ml[g]);
        if (bw_cov && bw_cov[g]) bmClose(bw_cov[g]);
        if (buffers) dm_pb_buf_free(&buffers[g]);
    }
    dm_sc_free_idmap(&idmap);
    free(idToGroup);
    free(outs);
    free(buffers);
    free(bw_ml);
    free(bw_cov);
    for (size_t i = 0; i < nMap; i++) free(mapping[i].cell);
    free(mapping);
    for (size_t i = 0; i < nGroups; i++) free(groups[i]);
    free(groups);
    free(input);
    free(outputPrefix);
    free(groupsPath);
    free(format);
    return 1;
}
