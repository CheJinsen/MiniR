/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000-2019 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
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
 *
 *  DESCRIPTION
 *
 *    The density of the lognormal distribution.
 */

#include "nmath.h"
#include "dpq.h"

double dlnorm(double x, double meanlog, double sdlog, bool give_log)
{
    if (ISNAN(x) || ISNAN(meanlog) || ISNAN(sdlog))
	   return x + meanlog + sdlog;

    if(sdlog < 0) ML_WARN_return_NAN;
    if(!R_FINITE(x) && log(x) == meanlog) return ML_NAN;/* log(x) - meanlog is NaN */
    if(sdlog == 0)
	   return (log(x) == meanlog) ? ML_POSINF : R_D__0;
    if(x <= 0) return R_D__0;

    double y = (log(x) - meanlog) / sdlog;
    return give_log ?
	    -(M_LN_SQRT_2PI   + 0.5 * y * y + log(x * sdlog)) :
	    M_1_SQRT_2PI * exp(-0.5 * y * y)  /	 (x * sdlog);
}
