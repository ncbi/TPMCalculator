#include <stdio.h>
#include <math.h>
#include "bmath.h"

double pt(double x, double n, int lower_tail, int log_p) {
    /* return  P[ T <= x ]	where
     * T ~ t_{n}  (t distrib. with n degrees of freedom).

     *	--> ./pnt.c for NON-central
     */
    double val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
#endif
    if (n <= 0.0) ML_ERR_return_NAN;

    if (!isfinite(x))
        return (x < 0) ? R_DT_0 : R_DT_1;
    if (!isfinite(n))
        return pnorm(x, 0.0, 1.0, lower_tail, log_p);

    nx = 1 + (x / n) * x;
    /* FIXME: This test is probably losing rather than gaining precision,
     * now that pbeta(*, log_p = TRUE) is much better.
     * Note however that a version of this test *is* needed for x*x > D_MAX */
    if (nx > 1e100) { /* <==>  x*x > 1e100 * n  */
        /* Danger of underflow. So use Abramowitz & Stegun 26.5.4
           pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
           with z = 1/nx,  a = n/2,  b= 1/2 :
         */
        double lval;
        lval = -0.5 * n * (2 * log(fabs(x)) - log(n))
                - lbeta(0.5 * n, 0.5) - log(0.5 * n);
        val = log_p ? lval : exp(lval);
    } else {
        val = (n > x * x)
                ? pbeta(x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
                : pbeta(1. / nx, n / 2., 0.5, /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
    if (x <= 0.)
        lower_tail = !lower_tail;

    if (log_p) {
        if (lower_tail) return log1p(-0.5 * exp(val));
        else return val - M_LN2; /* = log(.5* pbeta(....)) */
    } else {
        val /= 2.;
        return R_D_Cval(val);
    }
}
