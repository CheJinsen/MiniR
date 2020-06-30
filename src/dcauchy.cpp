/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Core Team
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
 *    The density of the Cauchy distribution.
 */

#include "nmath.h"
#include "dpq.h"

double dcauchy(double x, double location, double scale, bool give_log)
{
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(location) || ISNAN(scale))
		return x + location + scale;

    if (scale <= 0) ML_WARN_return_NAN;

    double y = (x - location) / scale;
    return give_log ?
		-log(M_PI * scale * (1. + y * y)) :
		1.0 / (M_PI * scale * (1. + y * y));
}