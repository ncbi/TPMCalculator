#include <stdio.h>
#include <math.h>
#include "bmath.h"

double pbeta_raw(double x, double a, double b, int lower_tail, int log_p) {
    // treat limit cases correctly here:
    if (a == 0 || b == 0 || !isfinite(a) || !isfinite(b)) {
        // NB:  0 < x < 1 :
        if (a == 0 && b == 0) // point mass 1/2 at each of {0,1} :
            return (log_p ? -M_LN2 : 0.5);
        if (a == 0 ) // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
            return R_DT_1;
        if (b == 0 ) // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
            return R_DT_0;
        // else, remaining case:  a = b = Inf : point mass 1 at 1/2
        if (x < 0.5) return R_DT_0;
        else return R_DT_1;
    }
    // Now:  0 < a < Inf;  0 < b < Inf

    double x1 = 0.5 - x + 0.5, w, wc;
    int ierr;
    //====
    bratio(a, b, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
    //====
    // ierr in {10,14} <==> bgrat() error code ierr-10 in 1:4; for 1 and 4, warned *there*
    if (ierr && ierr != 11 && ierr != 14)
        fprintf(stderr, "pbeta_raw(%g, a=%g, b=%g, ..) -> bratio() gave error code %d", x, a, b, ierr);
    return lower_tail ? w : wc;
} /* pbeta_raw() */

double pbeta(double x, double a, double b, int lower_tail, int log_p) {
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b)) return x + a + b;
#endif

    if (a < 0 || b < 0) ML_ERR_return_NAN;
    // allowing a==0 and b==0  <==> treat as one- or two-point mass

    if (x <= 0)
        return R_DT_0;
    if (x >= 1)
        return R_DT_1;

    return pbeta_raw(x, a, b, lower_tail, log_p);
}
