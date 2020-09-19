/*
 * This file is part of MiniR.
 * AUTHOR
 *   Catherine Loader, catherine@research.bell-labs.com.
 *   October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2015 The R Core Team
 * Copyright (C) 2005-2015 The R Foundation
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

double Binomial::cdf(double x, double n, double p, bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(n) || std::isnan(p))
		return x + n + p;
	if (!std::isfinite(n) || !std::isfinite(p)) return InfNaN::nan();

	bool nonint = fabs(n - nearbyint(n)) > 1e-7 * std::max(1.0, fabs(n));
	if (nonint) {
		std::cout << "non-integer n = " << n << std::endl;
		return InfNaN::nan();
	}
	n = nearbyint(n);
	if (n < 0 || p < 0 || p > 1) return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (x < 0) return R_DT_0;
	x = floor(x + 1e-7);
	if (n <= x) return R_DT_1;
	return Beta::cdf(p, x + 1, n - x, !lower_tail, log_p);
}

double Binomial::pdf(double x, double n, double p, bool log_p)
{
	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;

	/* NaNs propagated correctly */
	if (std::isnan(x) || std::isnan(n) || std::isnan(p)) return x + n + p;

	bool negInonint = n < 0.0 ||
		(fabs((n) - nearbyint(n)) > 1e-7 * std::max(1.0, fabs(n)));

	if (p < 0 || p > 1 || negInonint)
		return InfNaN::nan();

	bool nonint = fabs((x)-nearbyint(x)) > 1e-7 * std::max(1.0, fabs(x));
	if (nonint) {
		std::cout << "Warning: non-integer x = " << x << std::endl;
		return R_D__0;
	}

	if (x < 0 || !std::isfinite(x)) return R_D__0;

	n = nearbyint(n);
	x = nearbyint(x);

	return Base::binomialPdfRaw(x, n, p, 1 - p, log_p);
}

double Binomial::doSearch(double y, double* z, double p,
	double n, double pr, double incr)
{
	if (*z >= p) {

		for (;;) {
			double newz;
			if (y == 0 ||
				(newz = Binomial::cdf(y - incr, n, pr, true, false)) < p)
				return y;
			y = std::max(0.0, y - incr);
			*z = newz;
		}
	}
	else {
		for (;;) {
			y = std::min(y + incr, n);
			if (y == n ||
				(*z = Binomial::cdf(y, n, pr, true, false)) >= p)
				return y;
		}
	}
}

double Binomial::quantile(double p, double n, double pr,
	bool lower_tail, bool log_p)
{
	double q, mu, sigma, gamma, z, y;

	if (std::isnan(p) || std::isnan(n) || std::isnan(pr))
		return p + n + pr;

	if (!std::isfinite(n) || !std::isfinite(pr))
		return InfNaN::nan();
	/* if log_p is true, p = -Inf is a legitimate value */
	if (!std::isfinite(p) && !log_p)
		return InfNaN::nan();

	if (n != floor(n + 0.5)) return InfNaN::nan();
	if (pr < 0 || pr > 1 || n < 0)
		return InfNaN::nan();

	//R_Q_P01_boundaries(p, 0, n);
	if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? n : 0.0;
		if (p == InfNaN::neginf())
			return lower_tail ? 0.0 : n;
	}
	else {
		if (p < 0 || p > 1)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? 0.0 : n;
		if (p == 1)
			return lower_tail ? n : 0.0;
	}

	if (pr == 0.0 || n == 0.0) return 0.0;

	q = 1 - pr;
	if (q == 0.) return n;
	mu = n * pr;
	sigma = sqrt(n * pr * q);
	gamma = (q - pr) / sigma;

	double temp = log_p ? (lower_tail ? exp(p) : -expm1(p))
		: (lower_tail ? p : (0.5 - p + 0.5));
	if (!lower_tail || log_p) {
		p = temp;
		if (p == 0.0) return 0.0;
		if (p == 1.0) return n;
	}

	if (p + 1.01 * DBL_EPSILON >= 1.) return n;
	z = Normal::quantile(p, 0.0, 1.0, true, false);
	y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

	if (y > n)
		y = n;

	z = Binomial::cdf(y, n, pr, true, false);
	p *= 1 - 64 * DBL_EPSILON;

	if (n < 1e5) return doSearch(y, &z, p, n, pr, 1.0);
	{
		double incr = floor(n * 0.001), oldincr;
		do {
			oldincr = incr;
			y = doSearch(y, &z, p, n, pr, incr);
			incr = std::max(1.0, floor(incr / 100));
		} while (oldincr > 1 && incr > n * 1e-15);
		return y;
	}
}

int Binomial::rand(const int size, const double prob)
{
	if (!std::isfinite(size)) {
		return InfNaN::nan();
	}
	if (!std::isfinite(prob) || size < 0 || prob < 0.0 || prob > 1.0) {
		return InfNaN::nan();
	}

	std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::binomial_distribution<int> u(size, prob);
	return u(e);
}