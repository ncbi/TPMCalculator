#include <stdio.h>
#include <math.h>
#include "bmath.h"

/* These are recursive, so we should do a stack check */

#ifndef MATHLIB_STANDALONE
void R_CheckStack(void);
#endif

double lfastchoose(double n, double k) {
    return -1.0* log(n + 1.0) - lbeta(n - k + 1.0, k + 1.0);
}

/* mathematically the same:
   less stable typically, but useful if n-k+1 < 0 : */
static double lfastchoose2(double n, double k, int *s_choose) {
    double r;
    r = lgammafn_sign(n - k + 1., s_choose);
    return lgammafn(n + 1.) - lgammafn(k + 1.) - r;
}

double lchoose(double n, double k) {
    k = R_forceint(k);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n) || ISNAN(k)) return n + k;
#endif    
    if (k < 2) {
        if (k < 0) return -INFINITY;
        if (k == 0) return 0.;
        /* else: k == 1 */
        return log(fabs(n));
    }
    /* else: k >= 2 */
    if (n < 0) {
        return lchoose(-n + k - 1, k);
    } else if (R_IS_INT(n)) {
        n = R_forceint(n);
        if (n < k) return -INFINITY;
        /* k <= n :*/
        if (n - k < 2) return lchoose(n, n - k); /* <- Symmetry */
        /* else: n >= k+2 */
        return lfastchoose(n, k);
    }
    /* else non-integer n >= 0 : */
    if (n < k - 1) {
        int s;
        return lfastchoose2(n, k, &s);
    }
    return lfastchoose(n, k);
}

#define k_small_max 30

/* 30 is somewhat arbitrary: it is on the *safe* side:
 * both speed and precision are clearly improved for k < 30.
 */
double choose(double n, double k) {
    double r;    
    k = R_forceint(k);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n) || ISNAN(k)) return n + k;
#endif
    if (k < k_small_max) {
        int j;
        if (n - k < k && n >= 0 && R_IS_INT(n)) k = n - k; /* <- Symmetry */
        if (k < 0) return 0.;
        if (k == 0) return 1.;
        /* else: k >= 1 */
        r = n;
        for (j = 2; j <= k; j++)
            r *= (n - j + 1) / j;
        return R_IS_INT(n) ? R_forceint(r) : r;
        /* might have got rounding errors */
    }
    /* else: k >= k_small_max */
    if (n < 0) {
        r = choose(-n + k - 1, k);
        if (ODD(k)) r = -r;
        return r;
    } else if (R_IS_INT(n)) {
        n = R_forceint(n);        
        if (n < k) return 0.0;
        if (n - k < k_small_max) return choose(n, n - k); /* <- Symmetry */
        return R_forceint(exp(lfastchoose(n, k)));
    }
    /* else non-integer n >= 0 : */
    if (n < k - 1) {
        int s_choose;
        r = lfastchoose2(n, k, /* -> */ &s_choose);
        return s_choose * exp(r);
    }
    return exp(lfastchoose(n, k));
}
