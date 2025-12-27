#include "dmScShrinkage.h"
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

typedef struct {
    int row;
    int col;
    double val;
} dm_mtx_entry_t;

static int dm_mtx_entry_cmp(const void *a, const void *b) {
    const dm_mtx_entry_t *ea = a;
    const dm_mtx_entry_t *eb = b;
    if (ea->row != eb->row) return ea->row < eb->row ? -1 : 1;
    if (ea->col != eb->col) return ea->col < eb->col ? -1 : 1;
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

static void dm_sc_shrinkage_usage() {
    fprintf(stderr,
            "Usage: dmtools sc-shrinkage --mtx <mC.mtx> --cov <cov.mtx> --barcodes <barcodes.tsv> --features <features.tsv> --out <prefix> [options]\n"
            "Options:\n"
            "  --mtx <matrix.mtx>     MatrixMarket of methylated counts (mC)\n"
            "  --cov <cov.mtx>        MatrixMarket of coverage counts\n"
            "  --barcodes <tsv>       Cell barcodes\n"
            "  --features <tsv>       Region features\n"
            "  --out <prefix>         Output prefix\n"
            "  --prior-strength <M>   Prior strength M (alpha+beta) (default: 20)\n"
            "  --min-cov <N>          Minimum coverage per entry (default: 1)\n"
            "  --emit-var             Write posterior variance matrix\n"
            "  -h, --help             Show this help message\n");
}

int dm_sc_shrinkage_main(int argc, char **argv) {
    char *mtxPath = NULL;
    char *covPath = NULL;
    char *barcodesPath = NULL;
    char *featuresPath = NULL;
    char *outPrefix = NULL;
    double priorStrength = 20.0;
    int minCov = 1;
    int emitVar = 0;

    static struct option long_opts[] = {
        {"mtx", required_argument, 0, 1},
        {"cov", required_argument, 0, 2},
        {"barcodes", required_argument, 0, 3},
        {"features", required_argument, 0, 4},
        {"out", required_argument, 0, 5},
        {"prior-strength", required_argument, 0, 6},
        {"min-cov", required_argument, 0, 7},
        {"emit-var", no_argument, 0, 8},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int opt, long_idx = 0;
    while ((opt = getopt_long(argc, argv, "h", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 1: mtxPath = strdup(optarg); break;
            case 2: covPath = strdup(optarg); break;
            case 3: barcodesPath = strdup(optarg); break;
            case 4: featuresPath = strdup(optarg); break;
            case 5: outPrefix = strdup(optarg); break;
            case 6:
                priorStrength = atof(optarg);
                if (priorStrength <= 0.0) {
                    fprintf(stderr, "Error: --prior-strength must be positive.\n");
                    goto error;
                }
                break;
            case 7:
                minCov = atoi(optarg);
                if (minCov < 1) {
                    fprintf(stderr, "Error: --min-cov must be >= 1.\n");
                    goto error;
                }
                break;
            case 8: emitVar = 1; break;
            case 'h':
                dm_sc_shrinkage_usage();
                goto cleanup;
            default:
                dm_sc_shrinkage_usage();
                goto error;
        }
    }

    if (!mtxPath || !covPath || !barcodesPath || !featuresPath || !outPrefix) {
        fprintf(stderr, "Error: --mtx, --cov, --barcodes, --features, and --out are required.\n");
        dm_sc_shrinkage_usage();
        goto error;
    }

    dm_mtx_entry_t *mcEntries = NULL;
    dm_mtx_entry_t *covEntries = NULL;
    size_t nMc = 0, nCov = 0;
    int nRowsMc = 0, nColsMc = 0;
    int nRowsCov = 0, nColsCov = 0;
    if (dm_read_mtx(mtxPath, &nRowsMc, &nColsMc, &mcEntries, &nMc) != 0) goto error;
    if (dm_read_mtx(covPath, &nRowsCov, &nColsCov, &covEntries, &nCov) != 0) goto error;
    if (nRowsMc != nRowsCov || nColsMc != nColsCov) {
        fprintf(stderr, "Error: mC and cov matrices have different shapes.\n");
        goto error;
    }

    qsort(mcEntries, nMc, sizeof(dm_mtx_entry_t), dm_mtx_entry_cmp);
    qsort(covEntries, nCov, sizeof(dm_mtx_entry_t), dm_mtx_entry_cmp);

    double *sumCov = calloc(nColsMc, sizeof(double));
    double *sumMc = calloc(nColsMc, sizeof(double));
    if (!sumCov || !sumMc) goto error;

    size_t i = 0, j = 0;
    while (i < nMc && j < nCov) {
        int cmp = dm_mtx_entry_cmp(&mcEntries[i], &covEntries[j]);
        if (cmp == 0) {
            if (covEntries[j].val >= (double)minCov) {
                sumCov[mcEntries[i].col] += covEntries[j].val;
                sumMc[mcEntries[i].col] += mcEntries[i].val;
            }
            i++; j++;
        } else if (cmp < 0) {
            i++;
        } else {
            j++;
        }
    }

    double *priorAlpha = calloc(nColsMc, sizeof(double));
    double *priorBeta = calloc(nColsMc, sizeof(double));
    if (!priorAlpha || !priorBeta) goto error;
    for (int c = 0; c < nColsMc; c++) {
        if (sumCov[c] > 0.0) {
            double pbar = dm_sc_clamp_prob(sumMc[c] / sumCov[c]);
            priorAlpha[c] = pbar * priorStrength;
            priorBeta[c] = (1.0 - pbar) * priorStrength;
        }
    }

    size_t len = strlen(outPrefix) + 12;
    char *path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.mean.mtx", outPrefix);
    FILE *f = fopen(path, "w");
    free(path);
    if (!f) goto error;
    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d %d %zu\n", nRowsMc, nColsMc, nMc < nCov ? nMc : nCov);

    FILE *fvar = NULL;
    if (emitVar) {
        len = strlen(outPrefix) + 11;
        path = malloc(len);
        if (!path) goto error;
        snprintf(path, len, "%s.var.mtx", outPrefix);
        fvar = fopen(path, "w");
        free(path);
        if (!fvar) goto error;
        fprintf(fvar, "%%MatrixMarket matrix coordinate real general\n");
        fprintf(fvar, "%d %d %zu\n", nRowsMc, nColsMc, nMc < nCov ? nMc : nCov);
    }

    double *totalCov = calloc((size_t)nRowsMc, sizeof(double));
    double *totalMc = calloc((size_t)nRowsMc, sizeof(double));
    uint64_t *nFeat = calloc((size_t)nRowsMc, sizeof(uint64_t));
    if (!totalCov || !totalMc || !nFeat) goto error;

    i = 0;
    j = 0;
    while (i < nMc && j < nCov) {
        int cmp = dm_mtx_entry_cmp(&mcEntries[i], &covEntries[j]);
        if (cmp == 0) {
            if (covEntries[j].val >= (double)minCov) {
                double mean = dm_sc_posterior_mean(priorAlpha[mcEntries[i].col], priorBeta[mcEntries[i].col],
                                                   mcEntries[i].val, covEntries[j].val);
                fprintf(f, "%d %d %.6f\n", mcEntries[i].row + 1, mcEntries[i].col + 1, mean);
                if (fvar) {
                    double var = dm_sc_posterior_var(priorAlpha[mcEntries[i].col], priorBeta[mcEntries[i].col],
                                                     mcEntries[i].val, covEntries[j].val);
                    fprintf(fvar, "%d %d %.6f\n", mcEntries[i].row + 1, mcEntries[i].col + 1, var);
                }
                totalCov[mcEntries[i].row] += covEntries[j].val;
                totalMc[mcEntries[i].row] += mcEntries[i].val;
                nFeat[mcEntries[i].row] += 1;
            }
            i++; j++;
        } else if (cmp < 0) {
            i++;
        } else {
            j++;
        }
    }
    fclose(f);
    if (fvar) fclose(fvar);

    len = strlen(outPrefix) + 14;
    path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.barcodes.tsv", outPrefix);
    if (dm_sc_copy_file(barcodesPath, path) != 0) { free(path); goto error; }
    free(path);

    len = strlen(outPrefix) + 14;
    path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.features.tsv", outPrefix);
    if (dm_sc_copy_file(featuresPath, path) != 0) { free(path); goto error; }
    free(path);

    len = strlen(outPrefix) + 12;
    path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.qc.tsv", outPrefix);
    FILE *qc = fopen(path, "w");
    free(path);
    if (!qc) goto error;
    fprintf(qc, "cell_id\ttotal_cov\ttotal_mc\tn_features\n");
    FILE *barIn = fopen(barcodesPath, "r");
    if (!barIn) { fclose(qc); goto error; }
    char line[4096];
    size_t row = 0;
    while (fgets(line, sizeof(line), barIn)) {
        char *cell = strtok(line, "\t\n");
        if (!cell) continue;
        if (row >= (size_t)nRowsMc) break;
        fprintf(qc, "%s\t%.0f\t%.0f\t%" PRIu64 "\n",
                cell, totalCov[row], totalMc[row], nFeat[row]);
        row++;
    }
    fclose(barIn);
    fclose(qc);

    len = strlen(outPrefix) + 9;
    path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.mC.mtx", outPrefix);
    if (dm_sc_copy_file(mtxPath, path) != 0) { free(path); goto error; }
    free(path);

    len = strlen(outPrefix) + 10;
    path = malloc(len);
    if (!path) goto error;
    snprintf(path, len, "%s.cov.mtx", outPrefix);
    if (dm_sc_copy_file(covPath, path) != 0) { free(path); goto error; }
    free(path);

    free(mcEntries);
    free(covEntries);
    free(sumCov);
    free(sumMc);
    free(priorAlpha);
    free(priorBeta);
    free(totalCov);
    free(totalMc);
    free(nFeat);
    goto cleanup;

error:
    free(mcEntries);
    free(covEntries);
    free(sumCov);
    free(sumMc);
    free(priorAlpha);
    free(priorBeta);
    free(totalCov);
    free(totalMc);
    free(nFeat);
    return 1;

cleanup:
    free(mtxPath);
    free(covPath);
    free(barcodesPath);
    free(featuresPath);
    free(outPrefix);
    return 0;
}
