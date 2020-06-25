/*
 *  Mathlib : A C Library of Special Functions
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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cfloat>
#include <iostream>
#include <cstdbool>
#include <algorithm>

#include "Rmath.h"
//#include "R_ext/Random.h"

#define IEEE_754 1

#define R_forceint(x)   nearbyint(x)

//R >= 3.1.0: # define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7)
#define R_nonint(x) (fabs((x) - R_forceint(x)) > 1e-7*fmax2(1., fabs(x)))

#define _(String) (String)

#define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
#define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
#define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)
#define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) printf(fmt,x,x2,x3,x4,x5)

#define ISNAN(x) (std::isnan(x)!=0)

#define R_FINITE(x) std::isfinite(x)

#ifdef Windows
#define ML_POSINF	posinf()
#define ML_NEGINF	neginf()
#define ML_NAN		nan()

double f_div(const double a, const double b);
double posinf();
double neginf();
double nan();
#else
#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	(-1.0 / 0.0)
#define ML_NAN		(0.0 / 0.0)
#endif

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_WARN_return_NAN { ML_WARNING(ME_DOMAIN, ""); return ML_NAN; }

/* For a long time prior to R 2.3.0 ML_WARNING did nothing.
   We don't report ME_DOMAIN errors as the callers collect ML_NANs into
   a single warning.
 */
#define ML_WARNING(x, s) { \
   if(x > ME_DOMAIN) { \
       std::string msg = ""; \
       switch(x) { \
	   case ME_DOMAIN: \
		   msg = _("argument out of domain in ");	\
		   break; \
       case ME_RANGE: \
		   msg = _("value out of range in ");	\
		   break; \
       case ME_NOCONV: \
		   msg = _("convergence failed in ");	\
		   break; \
       case ME_PRECISION: \
		   msg = _("full precision may not have been achieved in "); \
		   break; \
       case ME_UNDERFLOW: \
		   msg = _("underflow occurred in ");	\
		   break; \
       } \
       std::cout << msg << s << std::endl; \
   } \
}

double pbeta_raw(double, double, double, bool, bool);

int	Rf_i1mach(int);

/* Used internally only */
double  Rf_d1mach(int);

void bratio(double a, double b, double x, double y,
	double* w, double* w1, int* ierr, bool log_p);

double bd0(double, double);
double stirlerr(double);  /* Stirling expansion "error" */
double lgammacor(double); /* log(gamma) correction */

	/* Chebyshev Series */
int chebyshev_init(double*, int, double);
double chebyshev_eval(double, const double*, const int);

double pgamma_raw(double, double, bool, bool);
