/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2018  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.

*/

#pragma once

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif
//#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
//# include <cmath>
//// See comment in R.h
//# ifdef __SUNPRO_CC
//using namespace std;
//# endif
//#else
//# include <math.h>
//#endif
//
//#ifdef NO_C_HEADERS
//# warning "use of NO_C_HEADERS is defunct and will be ignored"
//#endif
//
///*-- Mathlib as part of R --  define this for standalone : */
///* #undef MATHLIB_STANDALONE */
//
//#define R_VERSION_STRING "@PACKAGE_VERSION@"
//
//// Legacy defines -- C99 functions which R >= 3.5.0 reauires
//#ifndef HAVE_EXPM1
//# define HAVE_EXPM1 1
//#endif
//#ifndef HAVE_HYPOT
//# define HAVE_HYPOT 1
//#endif
//#ifndef HAVE_LOG1P
//# define HAVE_LOG1P 1
//#endif
//
//#ifndef HAVE_WORKING_LOG1P
//@RMATH_HAVE_WORKING_LOG1P@
//#endif
//
//#if !defined(HAVE_WORKING_LOG1P)
///* remap to avoid problems with getting the right entry point */
//double  Rlog1p(double);
//#define log1p Rlog1p
//#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


//#ifdef MATHLIB_STANDALONE
//# ifndef R_EXT_BOOLEAN_H_
///* "copy-paste" R_ext/Boolean.h if not already included: */
// #define R_EXT_BOOLEAN_H_
// #undef FALSE
// #undef TRUE
// typedef enum { FALSE = 0, TRUE } Rboolean;
//# endif
//#else
//# include <R_ext/Boolean.h>
//#endif
//
//
//#if !defined(MATHLIB_STANDALONE) && !defined(R_NO_REMAP_RMATH)
//#define bessel_i	Rf_bessel_i
//#define bessel_j	Rf_bessel_j
//#define bessel_k	Rf_bessel_k
//#define bessel_y	Rf_bessel_y
//#define bessel_i_ex	Rf_bessel_i_ex
//#define bessel_j_ex	Rf_bessel_j_ex
//#define bessel_k_ex	Rf_bessel_k_ex
//#define bessel_y_ex	Rf_bessel_y_ex
//#define beta		Rf_beta
//#define choose		Rf_choose
//#define dbeta		Rf_dbeta
//#define dbinom		Rf_dbinom
//#define dbinom_raw	Rf_dbinom_raw
//#define dcauchy		Rf_dcauchy
//#define dchisq		Rf_dchisq
//#define dexp		Rf_dexp
//#define df		Rf_df
//#define dgamma		Rf_dgamma
//#define dgeom		Rf_dgeom
//#define dhyper		Rf_dhyper
//#define digamma		Rf_digamma
//#define dlnorm		Rf_dlnorm
//#define dlogis		Rf_dlogis
//#define dnbeta		Rf_dnbeta
//#define dnbinom		Rf_dnbinom
//#define dnbinom_mu	Rf_dnbinom_mu
//#define dnchisq		Rf_dnchisq
//#define dnf		Rf_dnf
//#define dnorm4		Rf_dnorm4
//#define dnt		Rf_dnt
//#define dpois_raw	Rf_dpois_raw
//#define dpois		Rf_dpois
//#define dpsifn		Rf_dpsifn
//#define dsignrank	Rf_dsignrank
//#define dt		Rf_dt
//#define dtukey		Rf_dtukey
//#define dunif		Rf_dunif
//#define dweibull	Rf_dweibull
//#define dwilcox		Rf_dwilcox
//#define fmax2		Rf_fmax2
//#define fmin2		Rf_fmin2
//#define fprec		Rf_fprec
//#define fround		Rf_fround
//#define ftrunc		Rf_ftrunc
//#define fsign		Rf_fsign
//#define gammafn		Rf_gammafn
//#define imax2		Rf_imax2
//#define imin2		Rf_imin2
//#define lbeta		Rf_lbeta
//#define lchoose		Rf_lchoose
//#define lgammafn	Rf_lgammafn
//#define lgammafn_sign	Rf_lgammafn_sign
//#define lgamma1p	Rf_lgamma1p
//#define log1pexp       	Rf_log1pexp
//#define log1pmx		Rf_log1pmx
//#define logspace_add	Rf_logspace_add
//#define logspace_sub	Rf_logspace_sub
//#define logspace_sum	Rf_logspace_sum
//#define pbeta		Rf_pbeta
//#define pbeta_raw	Rf_pbeta_raw
//#define pbinom		Rf_pbinom
//#define pcauchy		Rf_pcauchy
//#define pchisq		Rf_pchisq
//#define pentagamma	Rf_pentagamma
//#define pexp		Rf_pexp
//#define pf		Rf_pf
//#define pgamma		Rf_pgamma
//#define pgeom		Rf_pgeom
//#define phyper		Rf_phyper
//#define plnorm		Rf_plnorm
//#define plogis		Rf_plogis
//#define pnbeta		Rf_pnbeta
//#define pnbinom		Rf_pnbinom
//#define pnbinom_mu     	Rf_pnbinom_mu
//#define pnchisq		Rf_pnchisq
//#define pnf		Rf_pnf
//#define pnorm5		Rf_pnorm5
//#define pnorm_both	Rf_pnorm_both
//#define pnt		Rf_pnt
//#define ppois		Rf_ppois
//#define psignrank	Rf_psignrank
//#define psigamma	Rf_psigamma
//#define pt		Rf_pt
//#define ptukey		Rf_ptukey
//#define punif		Rf_punif
//#define pweibull	Rf_pweibull
//#define pwilcox		Rf_pwilcox
//#define qbeta		Rf_qbeta
//#define qbinom		Rf_qbinom
//#define qcauchy		Rf_qcauchy
//#define qchisq		Rf_qchisq
//#define qchisq_appr	Rf_qchisq_appr
//#define qexp		Rf_qexp
//#define qf		Rf_qf
//#define qgamma		Rf_qgamma
//#define qgeom		Rf_qgeom
//#define qhyper		Rf_qhyper
//#define qlnorm		Rf_qlnorm
//#define qlogis		Rf_qlogis
//#define qnbeta		Rf_qnbeta
//#define qnbinom		Rf_qnbinom
//#define qnbinom_mu     	Rf_qnbinom_mu
//#define qnchisq		Rf_qnchisq
//#define qnf		Rf_qnf
//#define qnorm5		Rf_qnorm5
//#define qnt		Rf_qnt
//#define qpois		Rf_qpois
//#define qsignrank	Rf_qsignrank
//#define qt		Rf_qt
//#define qtukey		Rf_qtukey
//#define qunif		Rf_qunif
//#define qweibull	Rf_qweibull
//#define qwilcox		Rf_qwilcox
//#define rbeta		Rf_rbeta
//#define rbinom		Rf_rbinom
//#define rcauchy		Rf_rcauchy
//#define rchisq		Rf_rchisq
//#define rexp		Rf_rexp
//#define rf		Rf_rf
//#define rgamma		Rf_rgamma
//#define rgeom		Rf_rgeom
//#define rhyper		Rf_rhyper
//#define rlnorm		Rf_rlnorm
//#define rlogis		Rf_rlogis
//#define rmultinom	Rf_rmultinom
//#define rnbeta		Rf_rnbeta
//#define rnbinom		Rf_rnbinom
//#define rnbinom_mu     	Rf_rnbinom_mu
//#define rnchisq		Rf_rnchisq
//#define rnf		Rf_rnf
//#define rnorm		Rf_rnorm
//#define rnt		Rf_rnt
//#define rpois		Rf_rpois
//#define rsignrank	Rf_rsignrank
//#define rt		Rf_rt
//#define rtukey		Rf_rtukey
//#define runif		Rf_runif
//#define rweibull	Rf_rweibull
//#define rwilcox		Rf_rwilcox
//#define sign		Rf_sign
//#define tetragamma	Rf_tetragamma
//#define trigamma	Rf_trigamma
//#endif
//
//#define dnorm dnorm4
//#define pnorm pnorm5
//#define qnorm qnorm5
//
//#ifdef  __cplusplus
//extern "C" {
//#endif
//	/* R's versions with !R_FINITE checks */
//
double R_pow(double x, double y);
double R_pow_di(double, int);
//
//	/* Random Number Generators */
//
//double	norm_rand(void);
//double	unif_rand(void);
//double  R_unif_index(double);
//double	exp_rand(void);
//#ifdef MATHLIB_STANDALONE
//void	set_seed(unsigned int, unsigned int);
//void	get_seed(unsigned int *, unsigned int *);
//#endif
//
//	/* Normal Distribution */
//
double	dnorm(double x, double mean = 0.0, double sd = 1.0, bool log_p = false);
double	pnorm(double q, double mean = 0.0, double sd = 1.0,
	bool lower_tail = true, bool log_p = false);
double	qnorm(double p, double mean = 0.0, double sd = 1.0,
	bool lower_tail = true, bool log_p = false);
double	rnorm(double mu = 0.0, double sigma = 1.0);
void	pnorm_both(double, double *, double *, bool, bool);/* both tails */
//
//	/* Uniform Distribution */
//
double	dunif(double x, double min = 0.0, double max = 1.0, bool log = false);
double	punif(double q, double min = 0.0, double max = 1.0,
	bool lower_tail = true, bool log = false);
double	qunif(double p, double min = 0.0, double max = 1.0,
	bool lower_tail = true, bool log = false);
//double	runif(double, double);
//
	/* Gamma Distribution */
double	dgamma(double x, double shape, double scale, bool log);
double	pgamma(double q, double shape, double scale, bool lower_tail = true, bool log = false);
double	qgamma(double p, double shape, double scale, bool lower_tail = true, bool log = false);
//double	rgamma(double, double);
//
//double  log1pmx(double);
//double  log1pexp(double); // <-- ../nmath/plogis.c
double  lgamma1p(double);
double  logspace_add(double, double);
//double  logspace_sub(double, double);
//double  logspace_sum(const double *, int);
//
//	/* Beta Distribution */
//
double	dbeta(double x, double shape1, double shape2, bool log = false);
double	pbeta(double q, double shape1, double shape2,
	bool lower_tail = true, bool log_p = false);
double	qbeta(double p, double shape1, double shape2,
	bool lower_tail = true, bool log_p = false);
//double	rbeta(double, double);
//
//	/* Lognormal Distribution */
//
double	dlnorm(double x, double meanlog = 0.0, double sdlog = 1.0, bool log = false);
double	plnorm(double q, double meanlog = 0.0, double sdlog = 1.0,
	bool lower_tail = true, bool log = false);
double	qlnorm(double p, double meanlog = 0.0, double sdlog = 1.0,
	bool lower_tail = true, bool log = false);
//double	rlnorm(double, double);
//
//	/* Chi-squared Distribution */
//
double	dchisq(double x, double df, bool log = false);
double	pchisq(double q, double df, bool lower_tail = true, bool log = false);
double	qchisq(double p, double df, bool lower_tail = true, bool log = false);
//double	rchisq(double);
//
//	/* Non-central Chi-squared Distribution */
//
double	dnchisq(double, double, double, bool);
double	pnchisq(double, double, double, bool, bool);
double	qnchisq(double, double, double, bool, bool);
//double	rnchisq(double, double);
//
//	/* F Distibution */
//
double	df(double x, double df1, double df2, bool log = false);
double	pf(double q, double df1, double df2,
	bool lower_tail = true, bool log = false);
double	qf(double p, double df1, double df2,
	bool lower_tail = true, bool log = false);
//double	rf(double, double);
//
//	/* Student t Distibution */
//
double	dt(double x, double df, bool log);
double	pt(double q, double df, bool lower_tail = true, bool log_p = false);
double	qt(double p, double df, bool lower_tail = true, bool log_p = false);
//double	rt(double);
//
//	/* Binomial Distribution */
//
double  dbinom_raw(double x, double n, double p, double q, bool give_log);
double	dbinom(double x, double size, double prob, bool log = false);
double	pbinom(double q, double size, double prob,
	bool lower_tail = true, bool log_p = false);
double	qbinom(double p, double size, double prob,
	bool lower_tail = true, bool log_p = false);
//double	rbinom(double, double);
//
//	/* Multnomial Distribution */
//
//void	rmultinom(int, double*, int, int*);
//
//	/* Cauchy Distribution */
//
double	dcauchy(double x, double location = 0.0, double scale = 1.0, bool log = false);
double	pcauchy(double q, double location = 0.0, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
double	qcauchy(double p, double location = 0.0, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
//double	rcauchy(double, double);
//
//	/* Exponential Distribution */
//
double	dexp(double x, double rate = 1.0, bool log = false);
double	pexp(double q, double rate = 1.0,
	bool lower_tail = true, bool log_p = false);
double	qexp(double p, double rate = 1.0,
	bool lower_tail = true, bool log_p = false);
//double	rexp(double);
//
//	/* Geometric Distribution */
//
double	dgeom(double x, double prob, bool log = false);
double	pgeom(double q, double prob,
	bool lower_tail = true, bool log_p = false);
double	qgeom(double p, double prob,
	bool lower_tail = true, bool log_p = false);
//double	rgeom(double);
//
//	/* Hypergeometric Distibution */
//
double	dhyper(double x, double m, double n, double k, bool log = false);
double	phyper(double q, double m, double n, double k,
	bool lower_tail = true, bool log_p = false);
double	qhyper(double p, double m, double n, double k,
	bool lower_tail = true, bool log_p = false);
//double	rhyper(double, double, double);
//
//	/* Negative Binomial Distribution */
//
double	dnbinom(double x, double size, double prob, bool log = false);
double	pnbinom(double q, double size, double prob,
	bool lower_tail = true, bool log_p = false);
double	qnbinom(double p, double size, double prob,
	bool lower_tail = true, bool log_p = false);
//double	rnbinom(double, double);
//
double	dnbinom_mu(double, double, double, bool);
double	pnbinom_mu(double, double, double, bool, bool);
double	qnbinom_mu(double, double, double, bool, bool);
//double	rnbinom_mu(double, double);
//
//	/* Poisson Distribution */
//
double	dpois_raw (double x, double lambda, bool log);
double	dpois(double x, double lambda, bool log = false);
double	ppois(double q, double lambda, bool lower_tail = true, bool log = false);
double	qpois(double p, double lambda, bool lower_tail = true, bool log = false);
//double	rpois(double);
//
//	/* Weibull Distribution */
//
double	dweibull(double x, double shape, double scale = 1.0, bool log = false);
double	pweibull(double q, double shape, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
double	qweibull(double p, double shape, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
//double	rweibull(double, double);

//	/* Logistic Distribution */
//
double	dlogis(double x, double location = 0.0, double scale = 1.0, bool log = false);
double	plogis(double q, double location = 0.0, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
double	qlogis(double p, double location = 0.0, double scale = 1.0,
	bool lower_tail = true, bool log_p = false);
//double	rlogis(double, double);
//
//	/* Non-central Beta Distribution */
//
double	dnbeta(double, double, double, double, bool);
double	pnbeta(double, double, double, double, bool, bool);
double	qnbeta(double, double, double, double, bool, bool);
//double	rnbeta(double, double, double);
//
//	/* Non-central F Distribution */
//
double  dnf(double, double, double, double, bool);
double	pnf(double, double, double, double, bool, bool);
double	qnf(double, double, double, double, bool, bool);
//
//	/* Non-central Student t Distribution */
//
double	dnt(double, double, double, bool);
double	pnt(double, double, double, bool, bool);
double	qnt(double, double, double, bool, bool);
//
//	/* Studentized Range Distribution */
//
double	ptukey(double q, double nmeans, double df,
	double nranges = 1.0, bool lower_tail = true, bool log_p = false);
double	qtukey(double p, double nmeans, double df,
	double nranges = 1.0, bool lower_tail = true, bool log_p = false);
//
//	/* Wilcoxon Rank Sum Distribution */
//
double dwilcox(double x, double m, double n, bool log = false);
double pwilcox(double q, double m, double n,
	bool lower_tail = true, bool log_p = false);
double qwilcox(double p, double m, double n,
	bool lower_tail = true, bool log_p = false);
//double rwilcox(double, double);
//
	/* Wilcoxon Signed Rank Distribution */

double dsignrank(double x, double n, bool log = false);
double psignrank(double q, double n,
	bool lower_tail = true, bool log_p = false);
double qsignrank(double p, double n,
	bool lower_tail = true, bool log_p = false);
//double rsignrank(double);


/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);
//void    dpsifn(double, int, int, int, double*, int*, int*);
//double	psigamma(double, double);
//double	digamma(double);
//double	trigamma(double);
//double	tetragamma(double);
//double	pentagamma(double);
//
//double	beta(double, double);
double	lbeta(double, double);
//
double	choose(double, double);
double	lchoose(double, double);
//
//	/* Bessel Functions */
//
//double	bessel_i(double, double, double);
//double	bessel_j(double, double);
//double	bessel_k(double, double, double);
//double	bessel_y(double, double);
//double	bessel_i_ex(double, double, double, double *);
//double	bessel_j_ex(double, double, double *);
//double	bessel_k_ex(double, double, double, double *);
//double	bessel_y_ex(double, double, double *);
//
//
//	/* General Support Functions */
//
int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
//double	sign(double);
//double	fprec(double, double);
//double	fround(double, double);
//double	fsign(double, double);
//double	ftrunc(double);
//
//double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
//double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */
//
///* More accurate cos(pi*x), sin(pi*x), tan(pi*x)
//
//   These declarations might clash with system headers if someone had
//   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
//   defined (and we try, above).
//   We check for that via the value of __STDC_IEC_60559_FUNCS__
//*/
//#if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
double cospi(double);
double sinpi(double);
double tanpi(double);
//#endif
//
///* Compute the log of a sum or difference from logs of terms, i.e.,
// *
// *     log (exp (logx) + exp (logy))
// * or  log (exp (logx) - exp (logy))
// *
// * without causing overflows or throwing away too much accuracy:
// */
//double  logspace_add(double logx, double logy);
//double  logspace_sub(double logx, double logy);
//
//
///* ----------------- Private part of the header file ------------------- */
//
//#if defined(MATHLIB_STANDALONE) && !defined(MATHLIB_PRIVATE_H)
///* second is defined by nmath.h */
//
///* If isnan is a macro, as C99 specifies, the C++
//   math header will undefine it. This happens on macOS */
//# ifdef __cplusplus
//  int R_isnancpp(double); /* in mlutils.c */
//#  define ISNAN(x)     R_isnancpp(x)
//# else
//#  define ISNAN(x)     (isnan(x)!=0)
//# endif
//
//# define R_FINITE(x)    R_finite(x)
//int R_finite(double);
//
//# ifdef _WIN32  /* not Win32 as no config information */
//#  ifdef RMATH_DLL
//#   define R_EXTERN extern __declspec(dllimport)
//#  else
//#   define R_EXTERN extern
//#  endif
//R_EXTERN double NA_REAL;
//R_EXTERN double R_PosInf;
//R_EXTERN double R_NegInf;
//R_EXTERN int N01_kind;
//#  undef R_EXTERN
//#else
//extern int N01_kind;
//# endif
//
//#endif /* MATHLIB_STANDALONE */
//
//#ifdef  __cplusplus
//}
//#endif


