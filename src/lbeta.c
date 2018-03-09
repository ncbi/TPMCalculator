#include <stdio.h>
#include <math.h>
#include "bmath.h"

double lbeta(double a, double b) {
    double corr, p, q;

#ifdef IEEE_754
    if (ISNAN(a) || ISNAN(b))
        return a + b;
#endif
    p = q = a;
    if (b < p) p = b; /* := min(a,b) */
    if (b > q) q = b; /* := max(a,b) */

    /* both arguments must be >= 0 */
    if (p < 0)
        ML_ERR_return_NAN
    else if (p == 0) {
        return INFINITY;
    } else if (!isfinite(q)) { /* q == +Inf */
        return -INFINITY;
    }

    if (p >= 10) {
        /* p and q are big. */
        corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
        return log(q) * -0.5 + M_LN_SQRT_2PI + corr
                + (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    } else if (q >= 10) {
        /* p is small, but q is big. */
        corr = lgammacor(q) - lgammacor(p + q);
        double d = p + q;
        if (fabs(d) > 1e-25)
            d = log1p(-p / d);
        else
            return INFINITY;
        return lgammafn(p) + corr + p - p * log(p + q)
                + (q - 0.5) * d;
    } else
        /* p and q are small: p <= q < 10. */
        /* R change for very small args */
        if (p < 1e-306) return lgamma(p) + (lgamma(q) - lgamma(p + q));
    return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
}

