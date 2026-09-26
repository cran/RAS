/* Small dense OLS via Gauss-Jordan elimination with partial pivoting.
 * p is always tiny here (<=4: intercept/x/U/V or subsets), so plain O(p^3)
 * elimination on p<=8 is fine. Mirrors R's `solve(crossprod(X), crossprod(X,y))`
 * used throughout segmented's seg.lm.fit / davies.test internals. */
#ifndef RAS_SMALLLINALG_H
#define RAS_SMALLLINALG_H

#include <math.h>

/* In-place Gauss-Jordan on p x p row-major A, augmented with b (p) and,
 * if Ainv != NULL, the p x p identity. On success (0): b becomes the
 * solution x, and Ainv (if requested) becomes A^-1. Returns -1 if singular. */
static inline int ras_gauss_jordan_solve(double *A, double *b, int p, double *Ainv) {
    if (p > 8) return -1;
    double M[8][17];
    int width = p + 1 + (Ainv ? p : 0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) M[i][j] = A[i * p + j];
        M[i][p] = b[i];
        if (Ainv) for (int j = 0; j < p; j++) M[i][p + 1 + j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < p; col++) {
        int piv = col; double best = fabs(M[col][col]);
        for (int r = col + 1; r < p; r++) {
            double v = fabs(M[r][col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-12) return -1;
        if (piv != col) {
            for (int j = 0; j < width; j++) { double t = M[col][j]; M[col][j] = M[piv][j]; M[piv][j] = t; }
        }
        double d = M[col][col];
        for (int j = 0; j < width; j++) M[col][j] /= d;
        for (int r = 0; r < p; r++) {
            if (r == col) continue;
            double f = M[r][col];
            if (f == 0.0) continue;
            for (int j = 0; j < width; j++) M[r][j] -= f * M[col][j];
        }
    }
    for (int i = 0; i < p; i++) {
        b[i] = M[i][p];
        if (Ainv) for (int j = 0; j < p; j++) Ainv[i * p + j] = M[i][p + 1 + j];
    }
    return 0;
}

/* Unweighted OLS: X is n x p, column-major. coef[p] and *rss filled on
 * success (0). XtXinv (row-major p*p) filled too if non-NULL (needed for
 * davies.test's invXtX1[0,0]). Returns -1 if X'X is singular. */
static inline int ras_ols_fit(const double *X, const double *y, int n, int p,
                               double *coef, double *rss, double *XtXinv) {
    if (p > 8) return -1;
    double XtX[64], Xty[8];
    for (int i = 0; i < p; i++) {
        const double *xi = X + (size_t) i * n;
        double s = 0.0;
        for (int k = 0; k < n; k++) s += xi[k] * y[k];
        Xty[i] = s;
        for (int j = i; j < p; j++) {
            const double *xj = X + (size_t) j * n;
            double sij = 0.0;
            for (int k = 0; k < n; k++) sij += xi[k] * xj[k];
            XtX[i * p + j] = sij;
            XtX[j * p + i] = sij;
        }
    }
    double coef_local[8];
    for (int i = 0; i < p; i++) coef_local[i] = Xty[i];
    if (ras_gauss_jordan_solve(XtX, coef_local, p, XtXinv) != 0) return -1;
    double rss_ = 0.0;
    for (int k = 0; k < n; k++) {
        double fit = 0.0;
        for (int i = 0; i < p; i++) fit += X[(size_t) i * n + k] * coef_local[i];
        double r = y[k] - fit;
        rss_ += r * r;
    }
    for (int i = 0; i < p; i++) coef[i] = coef_local[i];
    *rss = rss_;
    return 0;
}

/* R's default (type-7) quantile algorithm, single probability, on an
 * ALREADY-ASCENDING-SORTED array xs of length n. */
static inline double ras_quantile7(const double *xs, int n, double p) {
    if (n <= 1) return xs[0];
    double h = (n - 1) * p; /* 0-based position */
    int lo = (int) floor(h);
    int hi = (int) ceil(h);
    if (lo < 0) lo = 0;
    if (hi > n - 1) hi = n - 1;
    if (lo == hi) return xs[lo];
    return xs[lo] + (h - lo) * (xs[hi] - xs[lo]);
}

static int ras_cmp_double(const void *a, const void *b) {
    double da = *(const double *) a, db = *(const double *) b;
    return (da > db) - (da < db);
}

#endif /* RAS_SMALLLINALG_H */
