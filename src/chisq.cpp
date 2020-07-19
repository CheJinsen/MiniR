/*
 * This file is part of MiniR.
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000 The R Core Team
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
 * along with MiniR. If not, see <https://www.gnu.org/licenses/>.
 */

#include "randist.h"
using namespace Randist;

double Chisq::pdf(double x, double df, bool give_log)
{
    return Gamma::pdf(x, df / 2.0, 2.0, give_log);
}

double Chisq::cdf(double x, double df, bool lower_tail, bool log_p)
{
    return Gamma::cdf(x, df / 2.0, 2.0, lower_tail, log_p);
}

double Chisq::quantile(double p, double df, bool lower_tail, bool log_p)
{
    return Gamma::quantile(p, 0.5 * df, 2.0, lower_tail, log_p);
}

double Chisq::rand(const double df)
{
	if (!std::isfinite(df) || df < 0.0) {
		return InfNaN::nan();
	}

	std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::chi_squared_distribution<double> u(df);
	return u(e);
}