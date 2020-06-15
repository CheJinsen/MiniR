/*
 *  Copyright (C) 1998-2020  The R Core Team
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
 *  https://www.R-project.org/Licenses/
 */

/* Private header file for use during compilation of Mathlib */

#pragma once

#include <cmath>
#include <iostream>
#include <cfloat>
#include <cstdbool>

int R_isnancpp(const double x);
#define ISNAN(x) R_isnancpp(x)

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)


#define ML_WARN_return_NAN { ML_WARNING(); return ML_NAN; }
void ML_WARNING();

#define R_FINITE(x) isfinite(x)	// isfinite defined in <cmath>
#define R_forceint(x) nearbyint(x)

/* Utilities for `dpq' handling (density/probability/quantile) */
/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define R_D_half (log_p ? -M_LN2 : 0.5)		// 1/2 (lower- or upper tail)


/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
// ==   R_DT_log when we already "know" log_p == TRUE


#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_WARN_return_NAN

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    ML_WARN_return_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    ML_WARN_return_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }

#define R_P_bounds_01(x, x_min, x_max)	\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1
/* is typically not quite optimal for (-Inf,Inf) where
 * you'd rather have */
#define R_P_bounds_Inf_01(x)			\
    if(!R_FINITE(x)) {				\
	if (x > 0) return R_DT_1;		\
	/* x < 0 */return R_DT_0;		\
    }



/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_nonint(x))

// for discrete d<distr>(x, ...) :
#define R_D_nonint_check(x)				\
   if(R_nonint(x)) {					\
       MATHLIB_WARNING(_("non-integer x = %f"), x);	\
	return R_D__0;					\
   }



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

