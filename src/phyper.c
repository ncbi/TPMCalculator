/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 1999-2012  The R Core Team
 *  Copyright (C) 2004	     Morten Welinder
 *  Copyright (C) 2004	     The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *	The distribution function of the hypergeometric distribution.
 *
 * Current implementation based on posting
 * From: Morten Welinder <terra@gnome.org>
 * Cc: R-bugs@biostat.ku.dk
 * Subject: [Rd] phyper accuracy and efficiency (PR#6772)
 * Date: Thu, 15 Apr 2004 18:06:37 +0200 (CEST)
 ......

 The current version has very serious cancellation issues.  For example,
 if you ask for a small right-tail you are likely to get total cancellation.
 For example,  phyper(59, 150, 150, 60, FALSE, FALSE) gives 6.372680161e-14.
 The right answer is dhyper(0, 150, 150, 60, FALSE) which is 5.111204798e-22.

 phyper is also really slow for large arguments.

 Therefore, I suggest using the code below. This is a sniplet from Gnumeric ...
 The code isn't perfect.  In fact, if  x*(NR+NB)  is close to	n*NR,
 then this code can take a while. Not longer than the old code, though.

 -- Thanks to Ian Smith for ideas.
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "bmath.h"

double dbinom_raw(double x, double n, double p, double q, int give_log) {
    double lf, lc;

    if (p == 0) return ((x == 0) ? R_D__1 : R_D__0);
    if (q == 0) return ((x == n) ? R_D__1 : R_D__0);

    if (x == 0) {
        if (n == 0) return R_D__1;
        lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q);
        return ( R_D_exp(lc));
    }
    if (x == n) {
        lc = (q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p);
        return ( R_D_exp(lc));
    }
    if (x < 0 || x > n) return ( R_D__0);

    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);

    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    lf = log(M_2PI) + log(x) + log1p(-x / n);

    return R_D_exp(lc - 0.5 * lf);
}

double dhyper(double x, double r, double b, double n, int give_log) {
    double p, q, p1, p2, p3;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(r) || ISNAN(b) || ISNAN(n))
        return x + r + b + n;
#endif

    if (R_D_negInonint(r) || R_D_negInonint(b) || R_D_negInonint(n) || n > r + b)
        ML_ERR_return_NAN;
    if (R_D_negInonint(x))
        return (0);

    x = R_D_forceint(x);
    r = R_D_forceint(r);
    b = R_D_forceint(b);
    n = R_D_forceint(n);

    if (n < x || r < x || n - x > b) return (R_D__0);
    if (n == 0) return ((x == 0) ? R_D__1 : R_D__0);

    p = ((double) n) / ((double) (r + b));
    q = ((double) (r + b - n)) / ((double) (r + b));

    p1 = dbinom_raw(x, r, p, q, give_log);
    p2 = dbinom_raw(n - x, b, p, q, give_log);
    p3 = dbinom_raw(n, r + b, p, q, give_log);

    return ( (give_log) ? p1 + p2 - p3 : p1 * p2 / p3);
}

static double pdhyper(double x, double NR, double NB, double n, int log_p) {
    /*
     * Calculate
     *
     *	    phyper (x, NR, NB, n, TRUE, FALSE)
     *   [log]  ----------------------------------
     *	       dhyper (x, NR, NB, n, FALSE)
     *
     * without actually calling phyper.  This assumes that
     *
     *     x * (NR + NB) <= n * NR
     *
     */
    double sum = 0;
    double term = 1;

    while (x > 0 && term >= DBL_EPSILON * sum) {
        term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
        sum += term;
        x--;
    }

    double ss = (double) sum;
    return log_p ? log1p(ss) : 1 + ss;
}

double phyper(double x, double NR, double NB, double n,
        int lower_tail, int log_p) {
    /* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */

    double d, pd;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(NR) || ISNAN(NB) || ISNAN(n))
        return x + NR + NB + n;
#endif

    x = floor(x + 1e-7);
    NR = R_D_forceint(NR);
    NB = R_D_forceint(NB);
    n = R_D_forceint(n);

    if (NR < 0 || NB < 0 || isinf(NR + NB) || n < 0 || n > NR + NB)
        ML_ERR_return_NAN;

    if (x * (NR + NB) > n * NR) {
        /* Swap tails.	*/
        double oldNB = NB;
        NB = NR;
        NR = oldNB;
        x = n - x - 1;
        lower_tail = !lower_tail;
    }

    if (x < 0)
        return R_DT_0;
    if (x >= NR || x >= n)
        return R_DT_1;

    d = dhyper(x, NR, NB, n, log_p);
    pd = pdhyper(x, NR, NB, n, log_p);

    return log_p ? R_DT_Log(d + pd) : R_D_Lval(d * pd);
}
