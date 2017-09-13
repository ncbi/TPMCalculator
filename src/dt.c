#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <limits.h>
#include "bmath.h"

double dt(double x, double n, int give_log) {
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
#endif
    if (n <= 0) ML_ERR_return_NAN;
    if (!R_FINITE(x))
        return R_D__0;
    if (!R_FINITE(n))
        return dnorm(x, 0., 1., give_log);

    double u, t = -bd0(n / 2., (n + 1) / 2.) + stirlerr((n + 1) / 2.) - stirlerr(n / 2.),
            x2n = x * x / n, // in  [0, Inf]
            ax = 0., // <- -Wpedantic
            l_x2n; // := log(sqrt(1 + x2n)) = log(1 + x2n)/2
    Rboolean lrg_x2n = (x2n > 1. / DBL_EPSILON);
    if (lrg_x2n) { // large x^2/n :
        ax = fabs(x);
        l_x2n = log(ax) - log(n) / 2.; // = log(x2n)/2 = 1/2 * log(x^2 / n)
        u = //  log(1 + x2n) * n/2 =  n * log(1 + x2n)/2 =
                n * l_x2n;
    } else if (x2n > 0.2) {
        l_x2n = log(1 + x2n) / 2.;
        u = n * l_x2n;
    } else {
        l_x2n = log1p(x2n) / 2.;
        u = -bd0(n / 2., (n + x * x) / 2.) + x * x / 2.;
    }

    //old: return  R_D_fexp(M_2PI*(1+x2n), t-u);

    // R_D_fexp(f,x) :=  (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
    // f = 2pi*(1+x2n)
    //  ==> 0.5*log(f) = log(2pi)/2 + log(1+x2n)/2 = log(2pi)/2 + l_x2n
    //	     1/sqrt(f) = 1/sqrt(2pi * (1+ x^2 / n))
    //		       = 1/sqrt(2pi)/(|x|/sqrt(n)*sqrt(1+1/x2n))
    //		       = M_1_SQRT_2PI * sqrt(n)/ (|x|*sqrt(1+1/x2n))
    if (give_log)
        return t - u - (M_LN_SQRT_2PI + l_x2n);

    // else :  if(lrg_x2n) : sqrt(1 + 1/x2n) ='= sqrt(1) = 1
    double I_sqrt_ = (lrg_x2n ? sqrt(n) / ax : exp(-l_x2n));
    return exp(t - u) * M_1_SQRT_2PI * I_sqrt_;
}
