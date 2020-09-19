/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2016 The R Core Team
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

double Poisson::pdf(double x, double lambda, bool give_log)
{
    if (std::isnan(x) || std::isnan(lambda))
        return x + lambda;

    if (lambda < 0) InfNaN::nan();
    //R_D_nonint_check(x);
	bool isnontint = fabs(x - nearbyint(x)) > 1e-7;
	if (isnontint) {
		std::cout << "Warning: non-integer x = " << x << std::endl;
		return give_log ? InfNaN::neginf() : 0.0;
	}

    if (x < 0 || !std::isfinite(x))
        return give_log ? InfNaN::neginf() : 0.0;

    x = nearbyint(x);
    return Base::poissonPdfRaw(x, lambda, give_log);
}

double Poisson::cdf(double x, double lambda, bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(lambda))
		return x + lambda;
	if (lambda < 0.0) return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (x < 0)		return R_DT_0;
	if (lambda == 0.)	return R_DT_1;
	if (!std::isfinite(x))	return R_DT_1;
	x = floor(x + 1e-7);

	return Gamma::cdf(lambda, x + 1, 1., !lower_tail, log_p);
}

double Poisson::doSearch(double y, double* z, double p,
	double lambda, double incr)
{
	if (*z >= p) {
		/* search to the left */
		for (;;) {
			if (y == 0 ||
				(*z = cdf(y - incr, lambda, true, false)) < p)
				return y;
			y = std::max(0.0, y - incr);
		}
	}
	else {		/* search to the right */

		for (;;) {
			y = y + incr;
			if ((*z = cdf(y, lambda, true, false)) >= p)
				return y;
		}
	}
}

double Poisson::quantile(double p, double lambda,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(lambda))
		return p + lambda;

	if (!std::isfinite(lambda))
		return InfNaN::nan();
	if (lambda < 0)
		return InfNaN::nan();
	//R_Q_P01_check(p);
	if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
		return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (lambda == 0) return 0;
	if (p == R_DT_0) return 0;
	if (p == R_DT_1) return InfNaN::posinf();

	double mu = lambda;
	double sigma = sqrt(lambda);
	double gamma = 1.0 / sigma;

	const double temp = lower_tail ? (p) : (0.5 - p +0.5);
	if (!lower_tail || log_p) {
		p = log_p ? (lower_tail ? exp(p) : -expm1(p)) : temp;
		if (p == 0.) return 0;
		if (p == 1.) return InfNaN::posinf();
	}

	if (p + 1.01 * DBL_EPSILON >= 1.) return InfNaN::posinf();

	double z = Normal::quantile(p, 0.0, 1.0, true, false);
	double y = nearbyint(mu + sigma * (z + gamma * (z * z - 1) / 6));
	z = cdf(y, lambda, true, false);

	/* fuzz to ensure left continuity; 1 - 1e-7 may lose too much : */
	p *= 1 - 64 * DBL_EPSILON;

	/* If the mean is not too large a simple search is OK */
	if (lambda < 1e5) return doSearch(y, &z, p, lambda, 1);
	/* Otherwise be a bit cleverer in the search */
	{
		double incr = floor(y * 0.001);
		double oldincr = 0.0;
		do {
			oldincr = incr;
			y = doSearch(y, &z, p, lambda, incr);
			incr = std::max(1.0, floor(incr / 100.0));
		} while (oldincr > 1 && incr > lambda * 1e-15);
		return y;
	}
}

int Poisson::rand(const double mu)
{
	if (!std::isfinite(mu) || mu < 0.0) {
		return InfNaN::nan();
	}

	std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::poisson_distribution<int> u(mu);
	return u(e);
}