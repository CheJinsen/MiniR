/*
 * This file is part of MiniR.
 * Copyright (C) 2020 Jinsen Che
 *
 * MiniR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MiniR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <cfloat>
#include "specfun.h"

namespace Randist
{
	class Base
	{
	protected:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double poissonPdfRaw(double x, double lambda, bool give_log);
		static double binomialPdfRaw(const double x, const double n, const double p,
			const double q, bool log_p);
		static long double nonCentralBetaCdfRaw(double x, double o_x, double a,
			double b, double ncp);
		static double nonCentralBetaCdf2(double x, double o_x, double a, double b, double ncp,
			bool lower_tail, bool log_p);
		static double tanpi(double x);
	};
}