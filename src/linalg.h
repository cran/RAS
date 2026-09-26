/* Small dense linear-algebra helpers shared by scan_c.c and gwas_c.c.
 * p (covariate count) is always small (~10-20), so plain loops are fine. */
#ifndef RAS_LINALG_H
#define RAS_LINALG_H

/* t[p] = Mmat(p x n, column-major) %*% g[n] */
static inline void ras_matvec_Mg(const double *Mmat, const double *g, int p, int n, double *t) {
    for (int i = 0; i < p; i++) {
        double s = 0.0;
        for (int k = 0; k < n; k++) s += Mmat[i + (size_t) k * p] * g[k];
        t[i] = s;
    }
}

/* z[n] = Zmat(n x p, column-major) %*% t[p] */
static inline void ras_matvec_Zt(const double *Zmat, const double *t, int n, int p, double *z) {
    for (int k = 0; k < n; k++) z[k] = 0.0;
    for (int i = 0; i < p; i++) {
        double ti = t[i];
        const double *col = Zmat + (size_t) i * n;
        for (int k = 0; k < n; k++) z[k] += col[k] * ti;
    }
}

#endif /* RAS_LINALG_H */
