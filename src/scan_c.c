#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>   /* R_Calloc / R_Free */
#include <Rmath.h>      /* pt(), pchisq() -- same underlying code R itself uses */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rasbin.h"
#include "linalg.h"

/* ── sub.seq / this.seq construction, mirroring R exactly ────────────────── */

/* Builds c(0, seq(min_ws, max_ws, by=skip2)) then appends max_ws if the last
 * element isn't already max_ws. Returns the count and fills *out (caller frees). */
static int build_sub_seq(int min_ws, int max_ws, int skip2, int **out) {
    int cap = 2 + (max_ws - min_ws) / (skip2 > 0 ? skip2 : 1) + 2;
    int *seq = (int *) R_Calloc((size_t) cap, int);
    int n = 0;
    seq[n++] = 0;
    for (int w = min_ws; w <= max_ws; w += skip2) seq[n++] = w;
    if (n == 0 || seq[n - 1] != max_ws) seq[n++] = max_ws;
    *out = seq;
    return n;
}

/* ── The fused chunked forward scan ───────────────────────────────────────
 *
 * Streams SNP-column chunks from a .rasbin file, builds the PGS contribution
 * for each column on the fly (geno * weight, NA -> 0, restricted to holdout
 * rows -- this is compute_pgs_matrix's job, done per-chunk so the full
 * pgs.mat is never materialized), and reproduces screen_forward_max_region's
 * incremental expanding-window accumulation and per-window test exactly.
 *
 * mode = "continuous": exact FWL simple regression on residuals.
 * mode = "score":      binary Rao score test against a cached null fit.
 * (The legacy per-window "glm" Wald path is not yet ported -- see plan.)
 */
SEXP RAS_fused_pgs_scan(
    SEXP path_sexp, SEXP weights_sexp, SEXP holdout_idx_sexp, SEXP mode_sexp,
    SEXP keep_mask_sexp, SEXP Zmat_sexp, SEXP Mmat_sexp, SEXP resid_sexp,
    SEXP w0_sexp, SEXP syy_sexp, SEXP df_sexp,
    SEXP skip1_sexp, SEXP skip2_sexp, SEXP min_ws_sexp, SEXP max_ws_sexp,
    SEXP chunk_snps_sexp
) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    const char *mode  = CHAR(STRING_ELT(mode_sexp, 0));
    int is_continuous = (strcmp(mode, "continuous") == 0);
    if (!is_continuous && strcmp(mode, "score") != 0)
        error("RAS_fused_pgs_scan: mode must be 'continuous' or 'score', got '%s'", mode);

    int64_t N          = XLENGTH(weights_sexp);
    int     n_holdout  = LENGTH(holdout_idx_sexp);
    const double *weights     = REAL(weights_sexp);
    const int    *holdout_idx = INTEGER(holdout_idx_sexp);   /* 1-based row idx into geno */
    const int    *keep_mask   = INTEGER(keep_mask_sexp);     /* length n_holdout, 0/1 */

    int p   = ncols(Zmat_sexp);
    int n_c = nrows(Zmat_sexp);
    const double *Zmat = REAL(Zmat_sexp);
    const double *Mmat = REAL(Mmat_sexp);
    const double *resid = REAL(resid_sexp);  /* ryC (continuous) or r0 (score), length n_c */
    const double *w0     = REAL(w0_sexp);     /* score only */
    double SyyC = asReal(syy_sexp);
    double dfC  = asReal(df_sexp);

    int skip1 = asInteger(skip1_sexp);
    int skip2 = asInteger(skip2_sexp);
    int min_ws = asInteger(min_ws_sexp);
    int max_ws = asInteger(max_ws_sexp);
    int chunk_snps = asInteger(chunk_snps_sexp);
    if (chunk_snps < 1) chunk_snps = 1;

    /* keep_pos: positions (0-based, into the n_holdout-length accumulator)
     * that survive complete.cases -- length must equal n_c. */
    int *keep_pos = (int *) R_Calloc((size_t) n_c, int);
    {
        int kk = 0;
        for (int i = 0; i < n_holdout; i++) {
            if (keep_mask[i]) {
                if (kk >= n_c) { R_Free(keep_pos); error("RAS_fused_pgs_scan: keep_mask has more TRUEs than nrow(Zmat)"); }
                keep_pos[kk++] = i;
            }
        }
        if (kk != n_c) { R_Free(keep_pos); error("RAS_fused_pgs_scan: keep_mask TRUE count (%d) != nrow(Zmat) (%d)", kk, n_c); }
    }

    int *sub_seq; int n_sub = build_sub_seq(min_ws, max_ws, skip2, &sub_seq);

    int n_grid = (int) ((N - 1) / skip1) + 1;
    SEXP profile_sexp = PROTECT(allocVector(REALSXP, n_grid));
    double *profile = REAL(profile_sexp);

    FILE *fp = fopen(path, "rb");
    if (!fp) { R_Free(keep_pos); R_Free(sub_seq); UNPROTECT(1); error("RAS_fused_pgs_scan: cannot open '%s'", path); }
    rasbin_header_t hdr;
    int rc = rasbin_read_header(fp, &hdr);
    if (rc != 0) { fclose(fp); R_Free(keep_pos); R_Free(sub_seq); UNPROTECT(1); error("RAS_fused_pgs_scan: bad .rasbin header"); }
    if (hdr.n_snps != N) { fclose(fp); R_Free(keep_pos); R_Free(sub_seq); UNPROTECT(1);
        error("RAS_fused_pgs_scan: weights length (%lld) != n_snps in file (%lld)",
              (long long) N, (long long) hdr.n_snps); }
    int64_t n_samples = hdr.n_samples;

    /* Sliding column-chunk cache: pgsbuf holds n_holdout x buf_ncols doubles,
     * covering absolute columns [buf_start, buf_end] (1-based, inclusive).
     * Reloaded (one seek+read) whenever a grid point's full window span
     * isn't already covered. */
    double *pgsbuf = NULL;
    int64_t buf_start = 0, buf_end = -1; /* empty */
    double *raw = NULL;
    int64_t raw_cap = 0;

    double *acc = (double *) R_Calloc((size_t) n_holdout, double);
    double *g   = (double *) R_Calloc((size_t) n_c, double);
    double *t   = (double *) R_Calloc((size_t) p, double);
    double *z   = (double *) R_Calloc((size_t) n_c, double);

    for (int gi = 0; gi < n_grid; gi++) {
        int64_t this_start = 1 + (int64_t) gi * skip1;

        int64_t lo = this_start - max_ws + 1; if (lo < 1) lo = 1;
        int64_t hi = this_start + max_ws - 1; if (hi > N) hi = N;

        int already_covered = (buf_end >= buf_start && lo >= buf_start && hi <= buf_end);
        if (!already_covered) {
            int64_t new_start = lo;
            int64_t new_end   = hi;
            if (new_end - new_start + 1 < chunk_snps) {
                new_end = new_start + chunk_snps - 1;
                if (new_end > N) new_end = N;
            }
            int64_t ncols_new = new_end - new_start + 1;

            if (ncols_new > raw_cap) {
                if (raw) R_Free(raw);
                raw = (double *) R_Calloc((size_t) (n_samples * ncols_new), double);
                raw_cap = ncols_new;
            }
            int rcr = rasbin_fread_columns(fp, &hdr, new_start, ncols_new, raw);
            if (rcr != 0) {
                fclose(fp);
                if (pgsbuf) R_Free(pgsbuf);
                if (raw) R_Free(raw);
                R_Free(acc); R_Free(g); R_Free(t); R_Free(z);
                R_Free(keep_pos); R_Free(sub_seq);
                UNPROTECT(1);
                error("RAS_fused_pgs_scan: chunk read failed for columns [%lld, %lld]",
                      (long long) new_start, (long long) new_end);
            }

            if (pgsbuf) R_Free(pgsbuf);
            pgsbuf = (double *) R_Calloc((size_t) (n_holdout * ncols_new), double);
            for (int64_t k = 0; k < ncols_new; k++) {
                double wgt = weights[new_start - 1 + k];
                const double *rawcol = raw + (size_t) k * n_samples;
                double *pgscol = pgsbuf + (size_t) k * n_holdout;
                for (int m = 0; m < n_holdout; m++) {
                    double v = rawcol[holdout_idx[m] - 1];
                    if (ISNAN(v)) v = 0.0;
                    pgscol[m] = v * wgt;
                }
            }
            buf_start = new_start;
            buf_end   = new_end;
        }

        /* ws = 0: single-column PGS */
        {
            const double *col0 = pgsbuf + (size_t) (this_start - buf_start) * n_holdout;
            memcpy(acc, col0, (size_t) n_holdout * sizeof(double));
        }
        int64_t left0 = this_start, right0 = this_start;

        double best_p = 1.0;
        for (int si = 0; si < n_sub; si++) {
            int ws = sub_seq[si];
            if (si > 0) {
                int64_t left1  = this_start - ws + 1; if (left1 < 1) left1 = 1;
                int64_t right1 = this_start + ws - 1; if (right1 > N) right1 = N;

                if (left1 != left0) {
                    for (int64_t j = left1; j <= left0 - 1; j++) {
                        const double *colj = pgsbuf + (size_t) (j - buf_start) * n_holdout;
                        for (int m = 0; m < n_holdout; m++) acc[m] += colj[m];
                    }
                }
                if (right1 != right0) {
                    for (int64_t j = right0 + 1; j <= right1; j++) {
                        const double *colj = pgsbuf + (size_t) (j - buf_start) * n_holdout;
                        for (int m = 0; m < n_holdout; m++) acc[m] += colj[m];
                    }
                }
                left0 = left1; right0 = right1;
            }

            for (int i = 0; i < n_c; i++) g[i] = acc[keep_pos[i]];
            ras_matvec_Mg(Mmat, g, p, n_c, t);
            ras_matvec_Zt(Zmat, t, n_c, p, z);

            double pval;
            if (is_continuous) {
                double Srxx = 0.0, Srxy = 0.0;
                for (int k = 0; k < n_c; k++) {
                    double rg = g[k] - z[k];
                    Srxx += rg * rg;
                    Srxy += rg * resid[k];
                }
                if (!R_FINITE(Srxx) || Srxx <= 0.0 || !(dfC > 0.0)) {
                    pval = 1.0;
                } else {
                    double beta = Srxy / Srxx;
                    double RSS  = SyyC - beta * Srxy; if (RSS < 0.0) RSS = 0.0;
                    double SE   = sqrt((RSS / dfC) / Srxx);
                    double tval = beta / SE;
                    pval = 2.0 * pt(-fabs(tval), dfC, /*lower_tail=*/1, /*log_p=*/0);
                }
            } else {
                double U = 0.0, V = 0.0;
                for (int k = 0; k < n_c; k++) {
                    double gperp = g[k] - z[k];
                    U += g[k] * resid[k];
                    V += w0[k] * gperp * gperp;
                }
                if (!R_FINITE(V) || V <= 0.0) {
                    pval = 1.0;
                } else {
                    double Tstat = (U * U) / V;
                    pval = pchisq(Tstat, 1.0, /*lower_tail=*/0, /*log_p=*/0);
                }
            }
            if (pval < best_p) best_p = pval;
        }

        profile[gi] = -log10(best_p);
    }

    fclose(fp);
    if (pgsbuf) R_Free(pgsbuf);
    if (raw) R_Free(raw);
    R_Free(acc); R_Free(g); R_Free(t); R_Free(z);
    R_Free(keep_pos); R_Free(sub_seq);

    UNPROTECT(1);
    return profile_sexp;
}
