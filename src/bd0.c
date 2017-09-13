#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <limits.h>
#include "bmath.h"

double bd0(double x, double np) {
    double ej, s, s1, v;
    int j;

    if (!R_FINITE(x) || !R_FINITE(np) || np == 0.0) ML_ERR_return_NAN;

    if (fabs(x - np) < 0.1 * (x + np)) {
        v = (x - np) / (x + np); // might underflow to 0
        s = (x - np) * v; /* s using v -- change by MM */
        if (fabs(s) < DBL_MIN) return s;
        ej = 2 * x*v;
        v = v*v;
        for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
            ej *= v; // = v^(2j+1)
            s1 = s + ej / ((j << 1) + 1);
            if (s1 == s) /* last term was effectively 0 */
                return s1;
            s = s1;
        }
    }
    /* else:  | x - np |  is not too small */
    return (x * log(x / np) + np - x);
}
