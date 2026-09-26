/* Brent's 1-D minimizer, ported from the classic Netlib fmin.f / Brent (1973)
 * algorithm -- the same algorithm R's own C-level `Brent_fmin` (called by
 * stats::optimize()) implements. Deterministic, no RNG, so this should track
 * R's optimize() output closely given an identical objective function. Used
 * to port segmented::seg.lm.fit's per-iteration `optimize(search.minOK, c(0,1), ...)`
 * step-length line search. */
#ifndef RAS_BRENT_H
#define RAS_BRENT_H

#include <math.h>
#include <float.h>

/* Minimizes f over [ax, bx] to the given tol. Mirrors R's Brent_fmin exactly
 * (golden-section bracketing + parabolic interpolation), so f is evaluated
 * with the same sequence of trial points R's optimize() would use. */
static inline double ras_brent_fmin(double ax, double bx,
                                     double (*f)(double, void *), void *info,
                                     double tol) {
    const double c = (3.0 - sqrt(5.0)) * 0.5; /* squared inverse golden ratio */
    double a = ax, b = bx;
    double eps = DBL_EPSILON;
    double tol1 = eps + 1.0;
    eps = sqrt(eps);

    double v = a + c * (b - a);
    double w = v, x = v;
    double d = 0.0, e = 0.0;
    double fx = f(x, info);
    double fv = fx, fw = fx;
    double tol3 = tol / 3.0;

    for (;;) {
        double xm = (a + b) * 0.5;
        tol1 = eps * fabs(x) + tol3;
        double t2 = tol1 * 2.0;
        if (fabs(x - xm) <= t2 - (b - a) * 0.5) break;

        double p = 0.0, q = 0.0, r = 0.0;
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.0;
            if (q > 0.0) p = -p; else q = -q;
            r = e; e = d;
        }

        double u;
        if (fabs(p) >= fabs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        } else {
            d = p / q;
            u = x + d;
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }

        if (fabs(d) >= tol1) u = x + d;
        else if (d > 0.0) u = x + tol1;
        else u = x - tol1;

        double fu = f(u, info);

        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw; w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return x;
}

#endif /* RAS_BRENT_H */
