#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>   /* R_Calloc / R_Free */
#include <Rmath.h>      /* pt(), pnorm(), gammafn() -- same code R itself uses */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "smalllinalg.h"
#include "brent.h"

/* ── Ported from segmented::davies.test / seg.lm.fit / seg.lm.fit.boot ──────
 *
 * Specialised (not general-purpose) port for the ONE call shape RAS actually
 * uses: get_break_points(x,y,t) always calls
 *   segmented(lm(y ~ x), seg.Z = ~x, npsi = 1)   [seg.control() defaults]
 *   davies.test(fit_lm)                          [all defaults]
 * Tracing segmented.lm's design-matrix construction for this exact shape
 * collapses it to XREG=[1,x] (2 cols), Z=[x] (1 col), unweighted/no-offset
 * (id.w.offs always TRUE -> only the unweighted seg.lm.fit branch is live).
 * See RAS_fast side-branch memory / plan file for the full derivation.
 */

/* ---- davies.test: shared per-candidate t-stat, then branch on n<=300 ---- */

static double ras_davies_test(const double *x, const double *y, int n) {
    if (n < 5) return 1.0; /* davies.test needs df.res=n-3>0 and k=10 candidates in (x[2],x[n-1]) */

    double *xs = (double *) R_Calloc((size_t) n, double);
    memcpy(xs, x, (size_t) n * sizeof(double));
    qsort(xs, (size_t) n, sizeof(double), ras_cmp_double);

    const int k = 10;
    double lo = xs[1], hi = xs[n - 2];
    R_Free(xs);
    if (!(hi > lo)) return 1.0; /* degenerate x range */

    double t_stat[10], rss_i[10];
    int valid[10];
    int df_res = n - 3;

    double *X3 = (double *) R_Calloc((size_t) n * 3, double); /* [X.psi, 1, x] col-major */
    for (int j = 0; j < n; j++) { X3[n + j] = 1.0; X3[2 * n + j] = x[j]; }

    for (int i = 0; i < k; i++) {
        double psi = lo + (hi - lo) * ((double) i / (double) (k - 1));
        for (int j = 0; j < n; j++) { double u = x[j] - psi; X3[j] = u > 0.0 ? u : 0.0; }
        double coef[3], rss, XtXinv[9];
        if (ras_ols_fit(X3, y, n, 3, coef, &rss, XtXinv) != 0 || !(rss > 0.0) || !R_FINITE(coef[0])) {
            valid[i] = 0; continue;
        }
        double denom = sqrt((rss / df_res) * XtXinv[0]);
        if (!(denom > 0.0) || !R_FINITE(denom)) { valid[i] = 0; continue; }
        t_stat[i] = coef[0] / denom;
        rss_i[i] = rss;
        valid[i] = 1;
    }
    R_Free(X3);

    /* filter to valid candidates, preserving order (matches R's RIS[!is.na(...)]) */
    double tv[10], rv[10]; int nv = 0;
    for (int i = 0; i < k; i++) if (valid[i]) { tv[nv] = t_stat[i]; rv[nv] = rss_i[i]; nv++; }
    if (nv == 0) return 1.0;

    double M = fabs(tv[0]);
    for (int i = 1; i < nv; i++) if (fabs(tv[i]) > M) M = fabs(tv[i]);

    double p;
    if (n <= 300) {
        double V = 0.0;
        double prev_asin = NA_REAL;
        for (int i = 0; i < nv; i++) {
            double Z = tv[i] * sqrt(rv[i] / df_res);
            double RIS2 = (Z * Z) / (Z * Z + rv[i]);
            if (RIS2 < 0.0) RIS2 = 0.0;
            if (RIS2 > 1.0) RIS2 = 1.0;
            double a = asin(sqrt(RIS2));
            if (i > 0) V += fabs(a - prev_asin);
            prev_asin = a;
        }
        double u = (M * M) / ((double) df_res + M * M);
        double approxx = V * (pow(1.0 - u, (df_res - 1) / 2.0) * gammafn(df_res / 2.0 + 0.5)) /
                          (2.0 * gammafn(df_res / 2.0) * sqrt(M_PI));
        double p_naiv = pt(M, (double) df_res, /*lower_tail=*/0, /*log_p=*/0);
        p = 2.0 * (p_naiv + approxx);
    } else {
        double V = 0.0;
        for (int i = 1; i < nv; i++) V += fabs(tv[i] - tv[i - 1]);
        double approxx = V * exp(-(M * M) / 2.0) / sqrt(8.0 * M_PI);
        double p_naiv = pnorm(M, 0.0, 1.0, /*lower_tail=*/0, /*log_p=*/0);
        p = 2.0 * (p_naiv + approxx);
    }
    if (!R_FINITE(p) || p < 0.0) p = 1.0;
    if (p > 1.0) p = 1.0;
    return p;
}

/* ---- slope_test: no-intercept OLS y ~ x - 1, one-tailed t-test ---- */

static double ras_slope_test(const double *x, const double *y, int n, int lower_tail) {
    double sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < n; i++) { sxx += x[i] * x[i]; sxy += x[i] * y[i]; }
    if (!(sxx > 0.0)) return 1.0;
    double beta = sxy / sxx;
    double rss = 0.0;
    for (int i = 0; i < n; i++) { double r = y[i] - beta * x[i]; rss += r * r; }
    int df = n - 1;
    if (df <= 0) return 1.0;
    double se = sqrt((rss / df) / sxx);
    if (!(se > 0.0)) return 1.0;
    double tval = beta / se;
    return pt(tval, (double) df, lower_tail, 0);
}

/* ---- Muggeo single-run fit (seg.lm.fit, unweighted/no-offset branch) ---- */

typedef struct {
    int success;
    double psi, rss;
    double intercept, slope_left, slope_right;
} seglm_t;

typedef struct { const double *x, *y; int n; double psi_old, psi_new; double *X3buf; } search_ctx_t;

static double ras_search_min_obj(double h, void *vctx) {
    search_ctx_t *ctx = (search_ctx_t *) vctx;
    double psi = ctx->psi_new * h + ctx->psi_old * (1.0 - h);
    int n = ctx->n;
    /* X3buf layout is [U(col0), intercept(col1, fixed), x(col2, fixed)] --
     * only col0 (U) varies per trial h, must NOT touch col2 (x). */
    for (int k = 0; k < n; k++) { double u = ctx->x[k] - psi; ctx->X3buf[k] = u > 0.0 ? u : 0.0; }
    double coef[3], rss;
    if (ras_ols_fit(ctx->X3buf, ctx->y, n, 3, coef, &rss, NULL) != 0) return 1e300;
    return rss;
}

static seglm_t ras_seg_lm_fit(const double *x, const double *y, int n, double psi0, double alpha) {
    const double h_step = 1.25, toll = 1e-5;
    const int it_max = 30;
    seglm_t res; memset(&res, 0, sizeof(res));

    double *xs = (double *) R_Calloc((size_t) n, double);
    memcpy(xs, x, (size_t) n * sizeof(double));
    qsort(xs, (size_t) n, sizeof(double), ras_cmp_double);
    double limZ_lo = ras_quantile7(xs, n, alpha);
    double limZ_hi = ras_quantile7(xs, n, 1.0 - alpha);
    R_Free(xs);

    double *X3 = (double *) R_Calloc((size_t) n * 3, double); /* [U,1,x] -- U in col0 to match davies layout convention (irrelevant, just consistent) */
    double *X4 = (double *) R_Calloc((size_t) n * 4, double); /* [1,x,U,V] */
    for (int k = 0; k < n; k++) {
        X3[n + k] = 1.0; X3[2 * n + k] = x[k];
        X4[k] = 1.0; X4[n + k] = x[k];
    }

    double psi = psi0;
    for (int k = 0; k < n; k++) { double u = x[k] - psi; X3[k] = u > 0.0 ? u : 0.0; }
    double coef3[3], rss3;
    if (ras_ols_fit(X3, y, n, 3, coef3, &rss3, NULL) != 0) goto cleanup;
    double L0 = rss3;
    double final_coef[3] = { coef3[1], coef3[2], coef3[0] }; /* reorder to intercept,x,U for clarity below */
    double final_psi = psi;

    search_ctx_t ctx; ctx.x = x; ctx.y = y; ctx.n = n; ctx.X3buf = X3;

    for (int it = 1; it <= it_max; it++) {
        for (int k = 0; k < n; k++) {
            double u = x[k] - psi; double uu = u > 0.0 ? u : 0.0;
            X4[2 * n + k] = uu;
            X4[3 * n + k] = (x[k] > psi) ? -1.0 : 0.0;
        }
        double coef4[4], rss4;
        if (ras_ols_fit(X4, y, n, 4, coef4, &rss4, NULL) != 0) goto cleanup;
        double beta_c = coef4[2], gamma_c = coef4[3];
        if (beta_c == 0.0 || gamma_c == 0.0 || !R_FINITE(beta_c) || !R_FINITE(gamma_c)) goto cleanup; /* fix.npsi=TRUE -> hard stop */

        double psi_old = psi;
        double psi_prop = psi_old + h_step * gamma_c / beta_c;
        if (psi_prop < limZ_lo) psi_prop = limZ_lo;
        if (psi_prop > limZ_hi) psi_prop = limZ_hi;

        double tol_it = 0.001 + (pow(DBL_EPSILON, 0.25) - 0.001) * ((double) (it - 1) / (double) (it_max - 1));

        ctx.psi_old = psi_old; ctx.psi_new = psi_prop;
        double h_opt = ras_brent_fmin(0.0, 1.0, ras_search_min_obj, &ctx, tol_it);
        double L1 = ras_search_min_obj(h_opt, &ctx);
        psi = psi_prop * h_opt + psi_old * (1.0 - h_opt);

        /* refit 3-col at final psi for this iteration's coefficients (X3 already holds U at psi via last obj eval == this psi) */
        double coef3b[3], rss3b;
        for (int k = 0; k < n; k++) { double u = x[k] - psi; X3[0 * n + k] = u > 0.0 ? u : 0.0; }
        if (ras_ols_fit(X3, y, n, 3, coef3b, &rss3b, NULL) != 0) goto cleanup;
        final_coef[0] = coef3b[1]; final_coef[1] = coef3b[2]; final_coef[2] = coef3b[0];
        final_psi = psi;

        double epsilon = (L0 - L1) / (fabs(L0) + 0.1);
        L0 = L1;
        if (fabs(epsilon) <= toll) break;
    }

    res.success = 1;
    res.psi = final_psi;
    res.rss = L0;
    res.intercept = final_coef[0];
    res.slope_left = final_coef[1];
    res.slope_right = final_coef[1] + final_coef[2];

cleanup:
    R_Free(X3); R_Free(X4);
    return res;
}

/* ---- bootstrap-restart wrapper (seg.lm.fit.boot) ---- */

static seglm_t ras_seg_lm_fit_boot(const double *x, const double *y, int n, double psi0, double alpha) {
    const int n_boot = 10, break_boot = 5;

    seglm_t incumbent = ras_seg_lm_fit(x, y, n, psi0, alpha);
    if (!incumbent.success) return incumbent;

    double *xb = (double *) R_Calloc((size_t) n, double);
    double *yb = (double *) R_Calloc((size_t) n, double);
    double recent[16]; int recent_n = 0;

    GetRNGstate();
    for (int kboot = 0; kboot < n_boot; kboot++) {
        for (int i = 0; i < n; i++) {
            /* Case-resample with replacement. This draws from R's own RNG
             * stream (via unif_rand(), respecting the caller's set.seed())
             * but deliberately does NOT replicate segmented::seg.lm.fit.boot's
             * exact seed-from-mean(y)-digits derivation or R's own sample()
             * rejection-sampling algorithm bit-for-bit -- user decision
             * 2026-07-21: keep the bootstrap-RESTART STRUCTURE (10 restarts,
             * case-resample, keep-best-RSS, early-stop after 5 non-improving
             * rounds) rather than matching R's specific draws.
             *
             * SPEED LEVER: this whole 10-restart loop is the single biggest
             * remaining cost in a per-window ras_detect call (each restart
             * re-runs the full Muggeo iteration twice -- once on the
             * resample, once refit on the full data). If ras_detect still
             * needs to be faster after this port (it's the current AoU
             * bottleneck relative to scan/gwas -- see project-aou-perf
             * memory), dropping n_boot here (down to a handful of restarts,
             * or to 0 for a single deterministic fit) is the place to look.
             * Left at the segmented::seg.control() default of 10 on purpose:
             * the user chose to preserve Muggeo's bootstrap-restart
             * robustness rather than trade it away by default. */
            int idx = (int) (unif_rand() * n);
            if (idx >= n) idx = n - 1;
            xb[i] = x[idx]; yb[i] = y[idx];
        }
        seglm_t oboot = ras_seg_lm_fit(xb, yb, n, incumbent.psi, alpha);
        double psi_start = oboot.success ? oboot.psi : incumbent.psi;
        seglm_t ofull = ras_seg_lm_fit(x, y, n, psi_start, alpha);
        if (ofull.success && ofull.rss <= incumbent.rss) incumbent = ofull;

        if (recent_n < 16) recent[recent_n++] = incumbent.rss;
        else { for (int j = 0; j < 15; j++) recent[j] = recent[j + 1]; recent[15] = incumbent.rss; }
        if (recent_n >= break_boot) {
            int flat = 1;
            for (int j = recent_n - break_boot + 1; j < recent_n; j++) {
                if (fabs(round((recent[j] - recent[j - 1]) * 1e6) / 1e6) > 0.0) { flat = 0; break; }
            }
            if (flat) break;
        }
    }
    PutRNGstate();

    R_Free(xb); R_Free(yb);
    return incumbent;
}

/* ---- get_break_points equivalent (degenerate guard + fit + davies test) ---- */

typedef struct { int has_break; int break_idx /*1-based*/; double p_value, slope_left, slope_right; } breakpt_t;

static breakpt_t ras_get_break_points(const double *x, const double *y, int t) {
    breakpt_t res; res.has_break = 0; res.break_idx = 0; res.p_value = 1.0; res.slope_left = 0.0; res.slope_right = 0.0;
    if (t < 5) return res;

    /* degenerate guard: resid_sigma from lm(y~x) vs sqrt(eps)*(mean(|y|)+1) */
    double *X2 = (double *) R_Calloc((size_t) t * 2, double);
    for (int i = 0; i < t; i++) { X2[i] = 1.0; X2[t + i] = x[i]; }
    double coef2[2], rss2;
    int ok2 = (ras_ols_fit(X2, y, t, 2, coef2, &rss2, NULL) == 0);
    R_Free(X2);
    if (!ok2) return res;
    double resid_sigma = (t > 2) ? sqrt(rss2 / (t - 2)) : NA_REAL;
    double mean_abs_y = 0.0;
    for (int i = 0; i < t; i++) mean_abs_y += fabs(y[i]);
    mean_abs_y /= t;
    double thresh = sqrt(DBL_EPSILON) * (mean_abs_y + 1.0);
    if (!R_FINITE(resid_sigma) || resid_sigma < thresh) return res;

    double alpha = fmax(0.05, 1.0 / t);
    double xmin = x[0], xmax = x[0];
    for (int i = 1; i < t; i++) { if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i]; }
    double psi0 = (xmin + xmax) / 2.0;

    seglm_t fit = ras_seg_lm_fit_boot(x, y, t, psi0, alpha);
    if (!fit.success) return res;

    int best_i = 0; double best_d = fabs(x[0] - fit.psi);
    for (int i = 1; i < t; i++) { double d = fabs(x[i] - fit.psi); if (d < best_d) { best_d = d; best_i = i; } }

    double p = ras_davies_test(x, y, t);

    res.has_break = 1;
    res.break_idx = best_i + 1;
    res.p_value = R_FINITE(p) ? p : 1.0;
    res.slope_left = fit.slope_left;
    res.slope_right = fit.slope_right;
    return res;
}

/* ---- get_local_maximum ---- */

static int ras_get_local_maximum(const double *y, int N, int x0_1based, int window_size) {
    int lower = x0_1based - window_size; if (lower < 1) lower = 1;
    int upper = x0_1based + window_size; if (upper > N) upper = N;
    int best_i = lower; double best_v = y[lower - 1];
    for (int i = lower; i <= upper; i++) if (y[i - 1] > best_v) { best_v = y[i - 1]; best_i = i; }
    return best_i;
}

/* ---- growable int/double vectors for the driver's output accumulation ---- */

typedef struct { double *data; int len, cap; } dvec_t;
static void dvec_push(dvec_t *v, double x) {
    if (v->len >= v->cap) {
        v->cap = v->cap ? v->cap * 2 : 64;
        v->data = (double *) R_Realloc(v->data, (size_t) v->cap, double);
    }
    v->data[v->len++] = x;
}

/* ── Main driver: ports ras_detect's while/for sliding-window loop exactly ── */

SEXP RAS_detect_fast(
    SEXP x_sexp, SEXP y_sexp, SEXP p_thresh_sexp, SEXP min_length_sexp,
    SEXP skip_sexp, SEXP window_size_sexp, SEXP slope_win_sexp,
    SEXP thresh_left_sexp, SEXP thresh_right_sexp
) {
    int N = LENGTH(x_sexp);
    const double *x = REAL(x_sexp);
    const double *y = REAL(y_sexp);
    double p_thresh = asReal(p_thresh_sexp);
    int min_length = asInteger(min_length_sexp);
    int skip = asInteger(skip_sexp);
    int window_size = asInteger(window_size_sexp);
    int slope_win = asInteger(slope_win_sexp);
    double thresh_left = asReal(thresh_left_sexp);
    double thresh_right = asReal(thresh_right_sexp);

    dvec_t tau_hats = {0}, p_values = {0}, slope_left = {0}, slope_right = {0}, slope_angle = {0};
    dvec_t all_cp = {0}, all_p = {0};

    int start_t = min_length - window_size;

    while (1) {
        int end_point = N - min_length + 1;
        if (start_t > end_point) break;

        /* this.seq <- seq(start_t, end_point, by=skip); append end_point if missing */
        int n_seq = (end_point - start_t) / skip + 1;
        if (n_seq < 1) n_seq = 1;
        int *seq = (int *) R_Calloc((size_t) (n_seq + 1), int);
        int m = 0;
        for (int v = start_t; v <= end_point; v += skip) seq[m++] = v;
        if (m == 0 || seq[m - 1] != end_point) seq[m++] = end_point;
        int seq_last = seq[m - 1];

        int accepted = 0;
        for (int si = 0; si < m; si++) {
            start_t = seq[si];
            int t = window_size;
            if (start_t < 0) { t = window_size + start_t; start_t = 1; }
            if (start_t > (N - window_size + 1)) { t = N - start_t + 1; }
            if (t < 5) continue;

            const double *xw = x + (start_t - 1);
            const double *yw = y + (start_t - 1);
            breakpt_t br = ras_get_break_points(xw, yw, t);

            int global_break = -1;
            if (br.has_break) {
                global_break = start_t - 1 + br.break_idx;
                dvec_push(&all_cp, (double) global_break);
                dvec_push(&all_p, br.p_value);
            } else {
                continue;
            }

            int v1 = br.slope_left > 0.0;
            int v2 = br.slope_right < 0.0;

            if (br.p_value <= p_thresh) {
                if (v1 && v2) {
                    int tau = global_break;
                    int lower = tau - slope_win; if (lower < 1) lower = 1;
                    int upper = tau + slope_win; if (upper > N) upper = N;

                    int n1 = tau - lower + 1;
                    double *xl = (double *) R_Calloc((size_t) n1, double);
                    double *yl = (double *) R_Calloc((size_t) n1, double);
                    for (int i = 0; i < n1; i++) { xl[i] = x[lower - 1 + i] - x[tau - 1]; yl[i] = y[lower - 1 + i] - y[tau - 1]; }
                    double p1 = ras_slope_test(xl, yl, n1, 0);
                    R_Free(xl); R_Free(yl);

                    int n2 = upper - tau + 1;
                    double *xr = (double *) R_Calloc((size_t) n2, double);
                    double *yr = (double *) R_Calloc((size_t) n2, double);
                    for (int i = 0; i < n2; i++) { xr[i] = x[tau - 1 + i] - x[tau - 1]; yr[i] = y[tau - 1 + i] - y[tau - 1]; }
                    double p2 = ras_slope_test(xr, yr, n2, 1);
                    R_Free(xr); R_Free(yr);

                    if (p1 < thresh_left && p2 < thresh_right) {
                        dvec_push(&tau_hats, (double) tau);
                        dvec_push(&p_values, br.p_value);
                        dvec_push(&slope_left, br.slope_left);
                        dvec_push(&slope_right, br.slope_right);
                        double ang1 = atan(br.slope_left) * 180.0 / M_PI;
                        double ang2 = atan(br.slope_right) * 180.0 / M_PI;
                        dvec_push(&slope_angle, ang2 - ang1 + 180.0);

                        start_t = tau;
                        accepted = 1;
                        break; /* exit inner for-loop, resume outer while from tau */
                    } else {
                        all_p.data[all_p.len - 1] = 1.0;
                    }
                } else {
                    all_p.data[all_p.len - 1] = 1.0;
                }
            }
        }
        R_Free(seq);
        if (!accepted && start_t == seq_last) break;
    }

    /* refine tau_hats to nearest local peak in y (window.size=50, matching get_local_maximum's default) */
    for (int i = 0; i < tau_hats.len; i++) {
        int refined = ras_get_local_maximum(y, N, (int) tau_hats.data[i], 50);
        tau_hats.data[i] = (double) refined;
    }

    const char *names[] = {"tau_hats", "p.values", "slope.left", "slope.right",
                             "all.changepoints", "all.p.values", "slope.angle", "previous_tau_hats", ""};
    SEXP out = PROTECT(mkNamed(VECSXP, names));

    SEXP r_tau = PROTECT(allocVector(REALSXP, tau_hats.len));
    if (tau_hats.len > 0) memcpy(REAL(r_tau), tau_hats.data, (size_t) tau_hats.len * sizeof(double));
    SET_VECTOR_ELT(out, 0, r_tau);

    SEXP r_pv = PROTECT(allocVector(REALSXP, p_values.len));
    if (p_values.len > 0) memcpy(REAL(r_pv), p_values.data, (size_t) p_values.len * sizeof(double));
    SET_VECTOR_ELT(out, 1, r_pv);

    SEXP r_sl = PROTECT(allocVector(REALSXP, slope_left.len));
    if (slope_left.len > 0) memcpy(REAL(r_sl), slope_left.data, (size_t) slope_left.len * sizeof(double));
    SET_VECTOR_ELT(out, 2, r_sl);

    SEXP r_sr = PROTECT(allocVector(REALSXP, slope_right.len));
    if (slope_right.len > 0) memcpy(REAL(r_sr), slope_right.data, (size_t) slope_right.len * sizeof(double));
    SET_VECTOR_ELT(out, 3, r_sr);

    SEXP r_acp = PROTECT(allocVector(REALSXP, all_cp.len));
    if (all_cp.len > 0) memcpy(REAL(r_acp), all_cp.data, (size_t) all_cp.len * sizeof(double));
    SET_VECTOR_ELT(out, 4, r_acp);

    SEXP r_ap = PROTECT(allocVector(REALSXP, all_p.len));
    if (all_p.len > 0) memcpy(REAL(r_ap), all_p.data, (size_t) all_p.len * sizeof(double));
    SET_VECTOR_ELT(out, 5, r_ap);

    SEXP r_ang = PROTECT(allocVector(REALSXP, slope_angle.len));
    if (slope_angle.len > 0) memcpy(REAL(r_ang), slope_angle.data, (size_t) slope_angle.len * sizeof(double));
    SET_VECTOR_ELT(out, 6, r_ang);

    SEXP r_prev = PROTECT(allocVector(REALSXP, tau_hats.len));
    if (tau_hats.len > 0) memcpy(REAL(r_prev), tau_hats.data, (size_t) tau_hats.len * sizeof(double));
    SET_VECTOR_ELT(out, 7, r_prev);

    R_Free(tau_hats.data); R_Free(p_values.data); R_Free(slope_left.data); R_Free(slope_right.data);
    R_Free(all_cp.data); R_Free(all_p.data); R_Free(slope_angle.data);

    UNPROTECT(9);
    return out;
}

/* ---- standalone exports for granular validation against R's davies.test / slope_test ---- */

SEXP RAS_davies_test(SEXP x_sexp, SEXP y_sexp) {
    int n = LENGTH(x_sexp);
    double p = ras_davies_test(REAL(x_sexp), REAL(y_sexp), n);
    return ScalarReal(p);
}

SEXP RAS_slope_test(SEXP x_sexp, SEXP y_sexp, SEXP lower_tail_sexp) {
    int n = LENGTH(x_sexp);
    double p = ras_slope_test(REAL(x_sexp), REAL(y_sexp), n, asLogical(lower_tail_sexp));
    return ScalarReal(p);
}

SEXP RAS_get_break_points(SEXP x_sexp, SEXP y_sexp, SEXP t_sexp) {
    int t = asInteger(t_sexp);
    breakpt_t br = ras_get_break_points(REAL(x_sexp), REAL(y_sexp), t);
    const char *names[] = {"break.points", "p.values", "slope.left", "slope.right", ""};
    SEXP out = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, br.has_break ? ScalarInteger(br.break_idx) : R_NilValue);
    SET_VECTOR_ELT(out, 1, ScalarReal(br.p_value));
    SET_VECTOR_ELT(out, 2, br.has_break ? ScalarReal(br.slope_left) : R_NilValue);
    SET_VECTOR_ELT(out, 3, br.has_break ? ScalarReal(br.slope_right) : R_NilValue);
    UNPROTECT(1);
    return out;
}

/* Validation-only: single deterministic Muggeo fit, no bootstrap restarts
 * (mirrors R's segmented(..., control=seg.control(n.boot=0)) for a controlled
 * comparison against the C port's core Newton+Brent iteration in isolation). */
SEXP RAS_seg_fit_single(SEXP x_sexp, SEXP y_sexp) {
    int n = LENGTH(x_sexp);
    const double *x = REAL(x_sexp), *y = REAL(y_sexp);
    double xmin = x[0], xmax = x[0];
    for (int i = 1; i < n; i++) { if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i]; }
    double alpha = fmax(0.05, 1.0 / n);
    seglm_t fit = ras_seg_lm_fit(x, y, n, (xmin + xmax) / 2.0, alpha);
    const char *names[] = {"success", "psi", "rss", "slope.left", "slope.right", ""};
    SEXP out = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, ScalarLogical(fit.success));
    SET_VECTOR_ELT(out, 1, ScalarReal(fit.psi));
    SET_VECTOR_ELT(out, 2, ScalarReal(fit.rss));
    SET_VECTOR_ELT(out, 3, ScalarReal(fit.slope_left));
    SET_VECTOR_ELT(out, 4, ScalarReal(fit.slope_right));
    UNPROTECT(1);
    return out;
}

SEXP RAS_seg_fit_boot(SEXP x_sexp, SEXP y_sexp, SEXP n_boot_ignored) {
    (void) n_boot_ignored;
    int n = LENGTH(x_sexp);
    const double *x = REAL(x_sexp), *y = REAL(y_sexp);
    double xmin = x[0], xmax = x[0];
    for (int i = 1; i < n; i++) { if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i]; }
    double alpha = fmax(0.05, 1.0 / n);
    seglm_t fit = ras_seg_lm_fit_boot(x, y, n, (xmin + xmax) / 2.0, alpha);
    const char *names[] = {"success", "psi", "rss", "slope.left", "slope.right", ""};
    SEXP out = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, ScalarLogical(fit.success));
    SET_VECTOR_ELT(out, 1, ScalarReal(fit.psi));
    SET_VECTOR_ELT(out, 2, ScalarReal(fit.rss));
    SET_VECTOR_ELT(out, 3, ScalarReal(fit.slope_left));
    SET_VECTOR_ELT(out, 4, ScalarReal(fit.slope_right));
    UNPROTECT(1);
    return out;
}
