#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bmath.h"

//static double ***w; /* to store  cwilcox(i,j,k) -> w[i][j][k] */
//static int allocated_m, allocated_n;

typedef struct wilcox_s{
    double ***w;
    int m,n;
}wilcox_t;

void w_free(wilcox_t *w) {
    int i, j;

    for (i = w->m; i >= 0; i--) {
        for (j = w->n; j >= 0; j--) {
            if (w->w[i][j] != 0)
                free((void *) w->w[i][j]);
        }
        free((void *) w->w[i]);
    }
    free((void *) w->w);
    w->w = NULL;
    w->m = w->n = 0;
}

void w_init_maybe(wilcox_t *wilcox, int m, int n) {
    int i, j;

    if (m > n) {
        i = n;
        n = m;
        m = i;
    }
    if (wilcox->w && (m > wilcox->m || n > wilcox->n))
        w_free(wilcox); /* zeroes w */

    if (!wilcox->w) { /* initialize w[][] */
        m = max(m, WILCOX_MAX);
        n = max(n, WILCOX_MAX);
        wilcox->w = (double ***) calloc((size_t) m + 1, sizeof (double **));
        for (i = 0; i <= m; i++) {
            wilcox->w[i] = (double **) calloc((size_t) n + 1, sizeof (double *));
            for(j = 0; j <= n; j++ ){
                wilcox->w[i][j] = 0;
            }
        }
        wilcox->m = m;
        wilcox->n = n;
    }
}

void w_free_maybe(wilcox_t *wilcox) {
    if (wilcox->m > WILCOX_MAX || wilcox->n > WILCOX_MAX)
        w_free(wilcox);
}

/* This counts the number of choices with statistic = k */
double cwilcox(wilcox_t *wilcox, int k, int m, int n) {
    int c, u, i, j, l;
    u = m * n;
    if (k < 0 || k > u)
        return (0);
    c = (int) (u / 2);
    if (k > c)
        k = u - k; /* hence  k <= floor(u / 2) */
    if (m < n) {
        i = m;
        j = n;
    } else {
        i = n;
        j = m;
    } /* hence  i <= j */

    if (j == 0) /* and hence i == 0 */
        return (k == 0);


    /* We can simplify things if k is small.  Consider the Mann-Whitney 
       definition, and sort y.  Then if the statistic is k, no more 
       than k of the y's can be <= any x[i], and since they are sorted 
       these can only be in the first k.  So the count is the same as
       if there were just k y's. 
     */
    if (j > 0 && k < j) return cwilcox(wilcox, k, i, k);

    if (wilcox->w[i][j] == 0) {
        wilcox->w[i][j] = (double *) calloc((size_t) c + 1, sizeof (double));
        for (l = 0; l <= c; l++)
            wilcox->w[i][j][l] = -1.0;
    }
    if (wilcox->w[i][j][k] < 0) {
        if (j == 0) /* and hence i == 0 */
            wilcox->w[i][j][k] = (k == 0);
        else {
            wilcox->w[i][j][k] = cwilcox(wilcox, k - j, i - 1, j) + cwilcox(wilcox, k, i, j - 1);
        }

    }
//    printf("cwilcox: %d %d %d %f\n", i, j, k, wilcox->w[i][j][k]);
    return (wilcox->w[i][j][k]);
}

double dwilcox(double x, double m, double n, int give_log) {
    double d;
    wilcox_t *wilcox = (wilcox_t *) malloc(sizeof(wilcox_t)); 
    wilcox->w = NULL;
    wilcox->m = 0;
    wilcox->n = 0;

    m = R_forceint(m);
    n = R_forceint(n);
    if (m <= 0 || n <= 0)
        ML_ERR_return_NAN;

    if (fabs(x - R_forceint(x)) > 1e-7)
        return (R_D__0);
    x = R_forceint(x);
    if ((x < 0) || (x > m * n))
        return (R_D__0);

    int mm = (int) m, nn = (int) n, xx = (int) x;
    w_init_maybe(wilcox, mm, nn);
    d = give_log ?
            log(cwilcox(wilcox, xx, mm, nn)) - lchoose(m + n, n) :
            cwilcox(wilcox, xx, mm, nn) / choose(m + n, n);
    w_free(wilcox); 
    free(wilcox);
    return (d);
}

/* args have the same meaning as R function pwilcox */
double pwilcox(double q, double m, double n, int lower_tail, int log_p) {
    int i;
    double c, p;
    wilcox_t *wilcox = (wilcox_t *) malloc(sizeof(wilcox_t)); 
    wilcox->w = NULL;
    wilcox->m = 0;
    wilcox->n = 0;

    if (!isfinite(m) || !isfinite(n))
        ML_ERR_return_NAN;
    m = R_forceint(m);
    n = R_forceint(n);
    if (m <= 0 || n <= 0)
        ML_ERR_return_NAN;

    q = floor(q + 1e-7);

    if (q < 0.0)
        return (R_DT_0);
    if (q >= m * n)
        return (R_DT_1);

    int mm = (int) m, nn = (int) n;
    w_init_maybe(wilcox, mm, nn);
    c = choose(m + n, n);
    p = 0.0;
    /* Use summation of probs over the shorter range */
    if (q <= (m * n / 2)) {
        for (i = 0; i <= q; i++) {
            p += (cwilcox(wilcox, i, mm, nn) / c);
        }
    } else {
        q = m * n - q;
        for (i = 0; i < q; i++)
            p += cwilcox(wilcox, i, mm, nn) / c;
        lower_tail = !lower_tail; /* p = 1 - p; */
    }
    w_free(wilcox); 
    free(wilcox);
    return p;
} /* pwilcox */

/* x is 'p' in R function qwilcox */

double qwilcox(double x, double m, double n, int lower_tail, int log_p) {
    double c, p;
    wilcox_t *wilcox = (wilcox_t *) malloc(sizeof(wilcox_t)); 
    wilcox->w = NULL;
    wilcox->m = 0;
    wilcox->n = 0;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(m) || ISNAN(n))
        return (x + m + n);
#endif
    if (!isfinite(x) || !isfinite(m) || !isfinite(n))
        ML_ERR_return_NAN;
    R_Q_P01_check(x);

    m = R_forceint(m);
    n = R_forceint(n);
    if (m <= 0 || n <= 0)
        ML_ERR_return_NAN;

    if (x == R_DT_0)
        return (0);
    if (x == R_DT_1)
        return (m * n);

    if (log_p || !lower_tail)
        x = R_DT_qIv(x); /* lower_tail,non-log "p" */

    int mm = (int) m, nn = (int) n;
    w_init_maybe(wilcox, mm, nn);
    c = choose(m + n, n);
    p = 0;
    int q = 0;
    if (x <= 0.5) {
        x = x - 10 * DBL_EPSILON;
        for (;;) {
            p += cwilcox(wilcox, q, mm, nn) / c;
            if (p >= x)
                break;
            q++;
        }
    } else {
        x = 1 - x + 10 * DBL_EPSILON;
        for (;;) {
            p += cwilcox(wilcox, q, mm, nn) / c;
            if (p > x) {
                q = (int) (m * n - q);
                break;
            }
            q++;
        }
    }
    w_free(wilcox); 
    free(wilcox);
    return (q);
}

double rwilcox(double m, double n) {
    int i, j, k, *x;
    double r;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(m) || ISNAN(n))
        return (m + n);
#endif
    m = R_forceint(m);
    n = R_forceint(n);
    if ((m < 0) || (n < 0))
        ML_ERR_return_NAN;

    if ((m == 0) || (n == 0))
        return (0);

    r = 0.0;
    k = (int) (m + n);
    x = (int *) calloc((size_t) k, sizeof (int));
#ifdef MATHLIB_STANDALONE
    if (!x) MATHLIB_ERROR(_("wilcox allocation error %d"), 4);
#endif
    for (i = 0; i < k; i++)
        x[i] = i;
    for (i = 0; i < n; i++) {
        j = (int) floor(k * unif_rand());
        r += x[j];
        x[j] = x[--k];
    }
    free(x);
    return (r - n * (n - 1) / 2);
}

