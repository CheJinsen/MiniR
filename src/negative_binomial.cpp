/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000 and Feb, 2001.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2016 The R Core Team
 * Copyright (C) 2005-2016 The R Foundation
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

double NegBinomial::pdf(double x, double size, double prob, bool give_log)
{
    if (std::isnan(x) || std::isnan(size) || std::isnan(prob))
        return x + size + prob;
    if (prob <= 0 || prob > 1 || size < 0)
        return InfNaN::nan();

    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
    double R_D__1 = give_log ? 0.0 : 1.0;

    bool isNonInt = fabs((x)-nearbyint(x)) > 1e-7 * std::max(1.0, fabs(x));
    if (isNonInt) {
        std::cout << "Warning: non-integer x = " << x << std::endl;
        return give_log ? InfNaN::neginf() : 0.0;
    }

    if (x < 0 || !std::isfinite(x)) return R_D__0;
    /* limiting case as size approaches zero is point mass at zero */
    if (x == 0 && size == 0) return R_D__1;
    x = nearbyint(x);
    if (!std::isfinite(size)) size = DBL_MAX;

    double ans = Base::binomialPdfRaw(size, x + size, prob, 1 - prob, give_log);
    double p = ((double)size) / (size + x);
    return((give_log) ? log(p) + ans : p * ans);
}

double NegBinomial::pdf_mu(double x, double size, double mu, bool give_log)
{
    if (std::isnan(x) || std::isnan(size) || std::isnan(mu)) {
        return x + size + mu;
    }
	if (mu < 0 || size < 0) {
		return InfNaN::nan();
	}

	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
	double R_D__1 = give_log ? 0.0 : 1.0;

	double isNonInt = fabs((x)-nearbyint(x)) >
		1e-7 * std::max(1.0, fabs(x));
	if (isNonInt) {
		std::cout << "Warning: non-integer x = " << x << std::endl;
		return R_D__0;
	}
    if (x < 0 || !std::isfinite(x)) return R_D__0;
    if (x == 0 && size == 0) return R_D__1;

    x = nearbyint(x);
    if (!std::isfinite(size)) { // limit case: Poisson
        return Base::poissonPdfRaw(x, mu, give_log);
    }

	if (x == 0) {
		double temp = size * (size < mu ? log(size / (size + mu))
			: log1p(-mu / (size + mu)));
		return give_log ? temp : exp(temp);
	}

    if (x < 1e-10 * size) {
        double p = (size < mu ? log(size / (1 + size / mu))
			: log(mu / (1 + mu / size)));
		double temp = x * p - mu - lgamma(x + 1) +
			log1p(x * (x - 1) / (2 * size));
        return give_log ? temp : exp(temp);
    }
    else {
		double p = ((double)size) / (size + x);
        double ans = Base::binomialPdfRaw(size, x + size, size / (size + mu),
			mu / (size + mu), give_log);
        return (give_log) ? log(p) + ans : p * ans;
    }
}

double NegBinomial::cdf(double x, double size, double prob,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(size) || std::isnan(prob))
		return x + size + prob;
	if (!std::isfinite(size) || !std::isfinite(prob))	return InfNaN::nan();
	if (size < 0 || prob <= 0 || prob > 1)	return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (size == 0)
		return (x >= 0) ? R_DT_1 : R_DT_0;

	if (x < 0) return R_DT_0;
	if (!std::isfinite(x)) return R_DT_1;
	x = floor(x + 1e-7);
	return Beta::cdf(prob, size, x + 1, lower_tail, log_p);
}

double NegBinomial::cdf_mu(double x, double size, double mu,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(size) || std::isnan(mu))
		return x + size + mu;
	if (!std::isfinite(mu))	return InfNaN::nan();
	if (size < 0 || mu < 0)	return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (size == 0)
		return (x >= 0) ? R_DT_1 : R_DT_0;

	if (x < 0) return R_DT_0;
	if (!std::isfinite(x)) return R_DT_1;
	if (!std::isfinite(size))
		return Poisson::cdf(x, mu, lower_tail, log_p);

	x = floor(x + 1e-7);
	/* return
	 * pbeta(pr, size, x + 1, lower_tail, log_p);  pr = size/(size + mu), 1-pr = mu/(size+mu)
	 *
	 *= pbeta_raw(pr, size, x + 1, lower_tail, log_p)
	 *            x.  pin   qin
	 *=  bratio (pin,  qin, x., 1-x., &w, &wc, &ierr, log_p),  and return w or wc ..
	 *=  bratio (size, x+1, pr, 1-pr, &w, &wc, &ierr, log_p) */
	{
		int ierr = 0;
		double w = 0.0, wc = 0.0;
		Toms::bratio(size, x + 1, size / (size + mu), mu / (size + mu), &w, &wc, &ierr, log_p);
		if (ierr)
			std::cout << "Warning: cdf_mu() -> bratio() gave error code " << ierr << std::endl;
		return lower_tail ? w : wc;
	}
}

double NegBinomial::doSearch(double y, double* z, double p, double n,
	double pr, double incr)
{
	if (*z >= p) {	/* search to the left */
		for (;;) {
			if (y == 0 ||
				(*z = cdf(y - incr, n, pr, true, false)) < p)
				return y;
			y = std::max(0.0, y - incr);
		}
	}
	else {		/* search to the right */
		for (;;) {
			y = y + incr;
			if ((*z = cdf(y, n, pr, true, false)) >= p)
				return y;
		}
	}
}

double NegBinomial::quantile(double p, double size, double prob,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(size) || std::isnan(prob))
		return p + size + prob;

	if (prob == 0 && size == 0) return 0;

	if (prob <= 0 || prob > 1 || size < 0) return InfNaN::nan();

	if (prob == 1 || size == 0) return 0;

	//R_Q_P01_boundaries(p, 0, ML_POSINF);
	if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? InfNaN::posinf() : 0.0;
		if (p == InfNaN::neginf())
			return lower_tail ? 0.0 : InfNaN::posinf();
	}
	else {
		if (p < 0 || p > 1)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? 0.0 : InfNaN::posinf();
		if (p == 1)
			return lower_tail ? InfNaN::posinf() : 0.0;
	}

	double Q = 1.0 / prob;
	double P = (1.0 - prob) * Q;
	double mu = size * P;
	double sigma = sqrt(size * P * Q);
	double gamma = (Q + P) / sigma;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	const double temp1 = lower_tail ? (p) : (0.5 - (p)+0.5);
	if (!lower_tail || log_p) {
		p = log_p ? (lower_tail ? exp(p) : -expm1(p)) : temp1;
		if (p == R_DT_0) return 0;
		if (p == R_DT_1) return InfNaN::posinf();
	}

	if (p + 1.01 * DBL_EPSILON >= 1.0) return InfNaN::posinf();

	double z = Normal::quantile(p, 0.0, 1.0, true, false);
	double y = nearbyint(mu + sigma * (z + gamma * (z * z - 1) / 6));

	z = cdf(y, size, prob, true, false);
	p *= 1 - 64 * DBL_EPSILON;

	/* If the C-F value is not too large a simple search is OK */
	if (y < 1e5)
		return doSearch(y, &z, p, size, prob, 1);

	{
		double incr = floor(y * 0.001), oldincr;
		do {
			oldincr = incr;
			y = doSearch(y, &z, p, size, prob, incr);
			incr = std::max(1.0, floor(incr / 100.0));
		} while (oldincr > 1 && incr > y * 1e-15);
		return y;
	}
}

double NegBinomial::quantile_mu(double p, double size, double mu,
	bool lower_tail, bool log_p)
{
	if (size == InfNaN::posinf())
		return(Poisson::quantile(p, mu, lower_tail, log_p));
	return quantile(p, size, size / (size + mu), lower_tail, log_p);
}

int NegBinomial::rand(const int size, const double prob)
{
    if(!std::isfinite(prob) || std::isnan(size) ||
    	size <= 0 || prob <= 0 || prob > 1) {
    	return InfNaN::nan();
    }

    std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::negative_binomial_distribution<int> u(size, prob);
	return u(e);
}