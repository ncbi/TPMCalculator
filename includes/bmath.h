/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bmath.h
 * Author: veraalva
 *
 * Created on February 18, 2016, 2:49 PM
 */

#ifndef BMATH_H
#define BMATH_H

#ifdef __cplusplus
extern "C" {
#endif

#define R_NaN NAN
#define R_PosInf INFINITY
#define R_NegInf -INFINITY

#define ML_POSINF R_PosInf
#define ML_NEGINF R_NegInf
#define ML_NAN NAN

#define ML_ERR_return_NAN { fprintf(stderr, "ML_ERR_return_NAN\n"); return NAN; }

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880168872420969808   /* sqrt(2)        */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32 5.656854249492380195206754896838 /* sqrt(32) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI 1.772453850905516027298167483341 /* sqrt(pi) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#ifndef M_LN2
#define M_LN2  0.693147180559945309417232121458 /* ln(2) */
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

#ifndef M_1_PI
#define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi           */
#endif
#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169163975144   /* pi/2           */
#endif

#ifndef M_2PI
#define M_2PI  6.283185307179586476925286766559 /* 2*pi */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2 0.225791352644727432363097614947 /* log(sqrt(pi/2)) */
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif

#ifndef M_LOG10_2
#define M_LOG10_2 0.301029995663981195213738894724 /* log10(2) */
#endif

#define give_log log_p
#define R_D__0 (log_p ? -INFINITY : 0.)  /* 0 */
#define R_D__1 (log_p ? 0. : 1.)   /* 1 */
#define R_DT_0 (lower_tail ? R_D__0 : R_D__1)  /* 0 */
#define R_DT_1 (lower_tail ? R_D__1 : R_D__0)  /* 1 */
#define R_D_nonint(x)    (fabs((x) - floor((x)+0.5)) > 1e-7)
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
#define R_DT_Log(p) (lower_tail? (p) : R_Log1_Exp(p))
#define R_D_Lval(p) (lower_tail ? (p) : (0.5 - (p) + 0.5)) /*  p  */
#define R_D_Cval(p) (lower_tail ? (0.5 - (p) + 0.5) : (p)) /*  1 - p */
#define R_DT_CIv(p) (log_p ? (lower_tail ? -expm1(p) : exp(p)) : R_D_Cval(p))
#define R_D_qIv(p) (log_p ? exp(p) : (p))  /*  p  in qF(p,..) */
#define R_D_exp(x) (log_p ?  (x)  : exp(x)) /* exp(x) */  
#define R_D_log(p) (log_p ?  (p)  : log(p)) /* log(p) */
#define R_Q_P01_check(p) if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)) ) ML_ERR_return_NAN
#define R_DT_qIv(p) (log_p ? (lower_tail ? exp(p) : - expm1(p)) : R_D_Lval(p))
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_forceint(x) floor((x) + 0.5)
#define R_nonint(x)    (fabs((x) - (int) x) > 1e-7)
#define R_IS_INT(x)  (!R_nonint(x))
#define ISNAN(x)     (isnan(x)!=0)
#define R_FINITE(x) isfinite(x)
#define Rboolean bool
#define ODD(_K_) ((_K_) != 2 * floor((_K_) / 2.))

    /* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50
#define k_small_max 30

#define max(a,b)  ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b)  ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

#define ML_VALID(x) (!ISNAN(x))

#define ME_NONE  0
    /*	no error */
#define ME_DOMAIN 1
    /*	argument out of domain */
#define ME_RANGE 2
    /*	value out of range */
#define ME_NOCONV 4
    /*	process did not converge */
#define ME_PRECISION 8
    /*	does not have "full" precision */
#define ME_UNDERFLOW 16
    /*	and underflow occured (important for IEEE)*/

#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  \
    if (log_p) {     \
        if(p > 0)     \
            ML_ERR_return_NAN;    \
        if(p == 0) /* upper bound*/   \
            return lower_tail ? _RIGHT_ : _LEFT_; \
        if(p == ML_NEGINF)    \
            return lower_tail ? _LEFT_ : _RIGHT_; \
    }       \
    else { /* !log_p */     \
        if(p < 0 || p > 1)    \
            ML_ERR_return_NAN;    \
        if(p == 0)     \
            return lower_tail ? _LEFT_ : _RIGHT_; \
        if(p == 1)     \
            return lower_tail ? _RIGHT_ : _LEFT_; \
    }

    extern void set_seed(unsigned int i1, unsigned int i2);

    extern void get_seed(unsigned int *i1, unsigned int *i2);

    extern double unif_rand(void);

    extern double lchoose(double n, double k);

    extern double choose(double n, double k);

    extern double stirlerr(double n);

    extern double lgammafn_sign(double x, int *sgn);

    extern double lgammafn(double x);

    extern double gammafn(double x);

    extern double lgammacor(double x);

    extern double chebyshev_eval(double x, const double *a, const int n);

    extern double phyper(double x, double NR, double NB, double n,
            int lower_tail, int log_p);

    extern double lbeta(double a, double b);

    extern double dwilcox(double x, double m, double n, int give_log);

    extern double pwilcox(double q, double m, double n, int lower_tail, int log_p);

    extern double qwilcox(double x, double m, double n, int lower_tail, int log_p);

    extern double rwilcox(double m, double n);

    extern double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);

    extern void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

    extern void bratio(double a, double b, double x, double y, double *w, double *w1,
            int *ierr, int log_p);

    extern double Rf_d1mach(int i);

    extern double pbeta(double x, double a, double b, int lower_tail, int log_p);

    extern double pt(double x, double n, int lower_tail, int log_p);

    extern double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);

    extern double dnorm4(double x, double mu, double sigma, int give_log);

    extern double dt(double x, double n, int give_log);

    extern double bd0(double x, double np);


#define pnorm pnorm5
#define qnorm qnorm5
#define dnorm dnorm4

#ifdef __cplusplus
}
#endif

#endif /* BMATH_H */

