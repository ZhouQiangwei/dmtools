#include "binaMeth.h"
#include "dmCommon.h"
#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct {
    char *id;
    uint64_t n_sites;
    uint64_t total_cov;
    double sum_meth_weighted;
} dm_sc_qc_entry_t;

static void dm_sc_qc_free_entries(dm_sc_qc_entry_t *entries, size_t n) {
    if (!entries) return;
    for (size_t i = 0; i < n; i++) {
        free(entries[i].id);
    }
    free(entries);
}

static int dm_sc_qc_find_or_add(dm_sc_qc_entry_t **entries, size_t *n, size_t *cap, const char *id) {
    for (size_t i = 0; i < *n; i++) {
        if (strcmp((*entries)[i].id, id) == 0) return (int)i;
    }
    if (*n == *cap) {
        size_t newCap = (*cap == 0) ? 16 : (*cap * 2);
        dm_sc_qc_entry_t *tmp = realloc(*entries, newCap * sizeof(dm_sc_qc_entry_t));
        if (!tmp) return -1;
        *entries = tmp;
        *cap = newCap;
    }
    (*entries)[*n].id = strdup(id);
    if (!(*entries)[*n].id) return -1;
    (*entries)[*n].n_sites = 0;
    (*entries)[*n].total_cov = 0;
    (*entries)[*n].sum_meth_weighted = 0.0;
    (*n)++;
    return (int)(*n - 1);
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
            "  cell_id  n_sites  total_coverage  mean_coverage  mean_meth\n");
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

    dm_sc_qc_entry_t *entries = NULL;
    size_t nEntries = 0, capEntries = 0;

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
                const char *id = o->entryid ? o->entryid[j] : NULL;
                if (!id) continue;
                int idx = dm_sc_qc_find_or_add(&entries, &nEntries, &capEntries, id);
                if (idx < 0) {
                    bmDestroyOverlappingIntervals(o);
                    bmClose(fp);
                    goto error_entries;
                }
                entries[idx].n_sites++;
                entries[idx].total_cov += o->coverage[j];
                entries[idx].sum_meth_weighted += ((double)o->value[j]) * o->coverage[j];
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
    fprintf(out, "cell_id\tn_sites\ttotal_coverage\tmean_coverage\tmean_meth\n");
    for (size_t i = 0; i < nEntries; i++) {
        double mean_cov = entries[i].n_sites ? (double)entries[i].total_cov / (double)entries[i].n_sites : 0.0;
        double mean_meth = entries[i].total_cov ? entries[i].sum_meth_weighted / (double)entries[i].total_cov : 0.0;
        fprintf(out, "%s\t%" PRIu64 "\t%" PRIu64 "\t%.6f\t%.6f\n",
                entries[i].id,
                entries[i].n_sites,
                entries[i].total_cov,
                mean_cov,
                mean_meth);
    }
    fclose(out);
    goto success;

error_entries:
    dm_sc_qc_free_entries(entries, nEntries);
error:
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPath);
cleanup:
    return 1;

success:
    dm_sc_qc_free_entries(entries, nEntries);
    for (size_t i = 0; i < nInputs; i++) free(inputs[i]);
    free(inputs);
    free(outputPath);
    return 0;
}

