#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>
#include <Rmath.h>
#include <math.h>
#include <string.h>
#include "rasbin.h"
#include "linalg.h"

/* ── Continuous branch: chunked, streams training-row values only ────────
 *
 * Mirrors compute_gwas_weights()'s continuous branch (R/gwas.R) exactly:
 * NA-free columns use the batched closed-form OLS; columns with an NA within
 * the training rows fall back to a per-column na.omit + refit, matching
 * lm() + na.omit() to the same precision the original R vectorization was
 * validated against.
 */
SEXP RAS_chunked_gwas_continuous(SEXP path_sexp, SEXP training_idx_sexp,
                                  SEXP y_sexp, SEXP chunk_snps_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    int m = LENGTH(training_idx_sexp);
    const int *training_idx = INTEGER(training_idx_sexp);
    const double *y = REAL(y_sexp);
    int chunk_snps = asInteger(chunk_snps_sexp);
    if (chunk_snps < 1) chunk_snps = 1;

    FILE *fp = fopen(path, "rb");
    if (!fp) error("RAS_chunked_gwas_continuous: cannot open '%s'", path);
    rasbin_header_t hdr;
    if (rasbin_read_header(fp, &hdr) != 0) { fclose(fp); error("RAS_chunked_gwas_continuous: bad .rasbin header"); }
    int64_t n_samples = hdr.n_samples, N = hdr.n_snps;

    /* Global y mean-centering, shared by every NA-free column. */
    double ybar = 0.0; for (int i = 0; i < m; i++) ybar += y[i]; ybar /= m;
    double *yc = (double *) R_Calloc((size_t) m, double);
    double Syy = 0.0;
    for (int i = 0; i < m; i++) { yc[i] = y[i] - ybar; Syy += yc[i] * yc[i]; }
    double df = (double) (m - 2);

    SEXP out = PROTECT(allocMatrix(REALSXP, (int) N, 4));
    double *coef = REAL(out);
    for (int64_t i = 0; i < N * 4; i++) coef[i] = NA_REAL;

    double *xi = (double *) R_Calloc((size_t) m, double);
    double *raw = NULL; int64_t raw_cap = 0;

    for (int64_t cs = 1; cs <= N; cs += chunk_snps) {
        int64_t ce = cs + chunk_snps - 1; if (ce > N) ce = N;
        int64_t ncols = ce - cs + 1;
        if (ncols > raw_cap) { if (raw) R_Free(raw); raw = (double *) R_Calloc((size_t) (n_samples * ncols), double); raw_cap = ncols; }
        if (rasbin_fread_columns(fp, &hdr, cs, ncols, raw) != 0) {
            fclose(fp); R_Free(yc); R_Free(xi); if (raw) R_Free(raw); UNPROTECT(1);
            error("RAS_chunked_gwas_continuous: chunk read failed [%lld,%lld]", (long long) cs, (long long) ce);
        }

        for (int64_t k = 0; k < ncols; k++) {
            int64_t j = cs + k;                 /* 1-based SNP index */
            const double *col = raw + (size_t) k * n_samples;
            int has_na = 0;
            for (int i = 0; i < m; i++) {
                double v = col[training_idx[i] - 1];
                xi[i] = v;
                if (ISNAN(v)) has_na = 1;
            }

            if (!has_na) {
                double xbar = 0.0; for (int i = 0; i < m; i++) xbar += xi[i]; xbar /= m;
                double Sxx = 0.0, Sxy = 0.0;
                for (int i = 0; i < m; i++) { double xc = xi[i] - xbar; Sxx += xc * xc; Sxy += xc * yc[i]; }
                if (R_FINITE(Sxx) && Sxx > 0.0 && df > 0.0) {
                    double beta = Sxy / Sxx;
                    double RSS = Syy - beta * Sxy; if (RSS < 0.0) RSS = 0.0;
                    double SE = sqrt((RSS / df) / Sxx);
                    double tval = beta / SE;
                    double *row = coef + (j - 1);
                    row[0]          = beta;
                    row[(size_t) N] = SE;
                    row[2 * N]      = tval;
                    row[3 * N]      = 2.0 * pt(-fabs(tval), df, 1, 0);
                }
            } else {
                int ni = 0;
                for (int i = 0; i < m; i++) if (!ISNAN(xi[i])) ni++;
                if (ni >= 3) {
                    double xbar = 0.0, ybar_k = 0.0;
                    for (int i = 0; i < m; i++) if (!ISNAN(xi[i])) { xbar += xi[i]; ybar_k += y[i]; }
                    xbar /= ni; ybar_k /= ni;
                    double sxx = 0.0, sxy = 0.0, syy = 0.0;
                    for (int i = 0; i < m; i++) {
                        if (ISNAN(xi[i])) continue;
                        double xk = xi[i] - xbar, yk = y[i] - ybar_k;
                        sxx += xk * xk; sxy += xk * yk; syy += yk * yk;
                    }
                    double dfi = (double) (ni - 2);
                    if (R_FINITE(sxx) && sxx > 0.0 && dfi > 0.0) {
                        double b = sxy / sxx;
                        double rss = syy - b * sxy; if (rss < 0.0) rss = 0.0;
                        double se = sqrt((rss / dfi) / sxx);
                        double tt = b / se;
                        double *row = coef + (j - 1);
                        row[0]          = b;
                        row[(size_t) N] = se;
                        row[2 * N]      = tt;
                        row[3 * N]      = 2.0 * pt(-fabs(tt), dfi, 1, 0);
                    }
                }
            }
        }
    }

    fclose(fp);
    R_Free(yc); R_Free(xi); if (raw) R_Free(raw);
    UNPROTECT(1);
    return out;
}

/* ── Binary branch, NA-free-within-training-subset fast path ─────────────
 *
 * Mirrors compute_gwas_weights()'s binary FWL "clean" block: residualises
 * each SNP column against the (precomputed once, in R) covariate projection
 * Mmat = (W'W)^-1 W', then a simple regression on residuals. Columns with an
 * NA among train_ok_idx rows are left NA here and returned in `na_cols` for
 * an R-side fallback identical to the original per-column na.omit path
 * (which needs model.matrix()'s factor handling and so isn't ported to C).
 */
SEXP RAS_chunked_gwas_binary_clean(SEXP path_sexp, SEXP train_ok_idx_sexp,
                                    SEXP Mmat_sexp, SEXP Wmat_sexp,
                                    SEXP ryb_sexp, SEXP Syy_sexp, SEXP df_sexp,
                                    SEXP chunk_snps_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    int nB = LENGTH(train_ok_idx_sexp);
    const int *train_ok_idx = INTEGER(train_ok_idx_sexp);
    int p = ncols(Wmat_sexp);
    const double *Mmat = REAL(Mmat_sexp);
    const double *Wmat = REAL(Wmat_sexp);
    const double *ryb  = REAL(ryb_sexp);
    double Syy = asReal(Syy_sexp);
    double dfree = asReal(df_sexp);
    int chunk_snps = asInteger(chunk_snps_sexp);
    if (chunk_snps < 1) chunk_snps = 1;

    FILE *fp = fopen(path, "rb");
    if (!fp) error("RAS_chunked_gwas_binary_clean: cannot open '%s'", path);
    rasbin_header_t hdr;
    if (rasbin_read_header(fp, &hdr) != 0) { fclose(fp); error("RAS_chunked_gwas_binary_clean: bad .rasbin header"); }
    int64_t n_samples = hdr.n_samples, N = hdr.n_snps;

    SEXP coef_sexp = PROTECT(allocMatrix(REALSXP, (int) N, 4));
    double *coef = REAL(coef_sexp);
    for (int64_t i = 0; i < N * 4; i++) coef[i] = NA_REAL;

    /* Overallocate; trim to the actual count before returning. */
    int *na_cols_buf = (int *) R_Calloc((size_t) N, int);
    int n_na_cols = 0;

    double *g  = (double *) R_Calloc((size_t) nB, double);
    double *t  = (double *) R_Calloc((size_t) p, double);
    double *z  = (double *) R_Calloc((size_t) nB, double);
    double *raw = NULL; int64_t raw_cap = 0;

    for (int64_t cs = 1; cs <= N; cs += chunk_snps) {
        int64_t ce = cs + chunk_snps - 1; if (ce > N) ce = N;
        int64_t ncols_chunk = ce - cs + 1;
        if (ncols_chunk > raw_cap) { if (raw) R_Free(raw); raw = (double *) R_Calloc((size_t) (n_samples * ncols_chunk), double); raw_cap = ncols_chunk; }
        if (rasbin_fread_columns(fp, &hdr, cs, ncols_chunk, raw) != 0) {
            fclose(fp); R_Free(g); R_Free(t); R_Free(z); R_Free(na_cols_buf); if (raw) R_Free(raw);
            UNPROTECT(1);
            error("RAS_chunked_gwas_binary_clean: chunk read failed [%lld,%lld]", (long long) cs, (long long) ce);
        }

        for (int64_t k = 0; k < ncols_chunk; k++) {
            int64_t j = cs + k;
            const double *col = raw + (size_t) k * n_samples;
            int has_na = 0;
            for (int i = 0; i < nB; i++) {
                double v = col[train_ok_idx[i] - 1];
                g[i] = v;
                if (ISNAN(v)) { has_na = 1; break; }
            }
            if (has_na) {
                na_cols_buf[n_na_cols++] = (int) j;
                continue;
            }

            ras_matvec_Mg(Mmat, g, p, nB, t);
            ras_matvec_Zt(Wmat, t, nB, p, z);
            double Srxx = 0.0, Srxy = 0.0;
            for (int i = 0; i < nB; i++) {
                double rg = g[i] - z[i];
                Srxx += rg * rg;
                Srxy += rg * ryb[i];
            }
            if (R_FINITE(Srxx) && Srxx > 0.0 && dfree > 0.0) {
                double beta = Srxy / Srxx;
                double RSS = Syy - beta * Srxy; if (RSS < 0.0) RSS = 0.0;
                double SE = sqrt((RSS / dfree) / Srxx);
                double tval = beta / SE;
                double *row = coef + (j - 1);
                row[0]          = beta;
                row[(size_t) N] = SE;
                row[2 * N]      = tval;
                row[3 * N]      = 2.0 * pt(-fabs(tval), dfree, 1, 0);
            }
        }
    }

    fclose(fp);
    R_Free(g); R_Free(t); R_Free(z); if (raw) R_Free(raw);

    SEXP na_cols_sexp = PROTECT(allocVector(INTSXP, n_na_cols));
    if (n_na_cols > 0) memcpy(INTEGER(na_cols_sexp), na_cols_buf, (size_t) n_na_cols * sizeof(int));
    R_Free(na_cols_buf);

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, coef_sexp);
    SET_VECTOR_ELT(result, 1, na_cols_sexp);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("coef_mat"));
    SET_STRING_ELT(names, 1, mkChar("na_cols"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(4);
    return result;
}
