/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2015 The R Core Team
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

double NonCentralBeta::pdf(double x, double a, double b, double ncp, bool give_log)
{
    constexpr double eps = 1.e-15;

    if (std::isnan(x) || std::isnan(a) || std::isnan(b) || std::isnan(ncp)) {
        return x + a + b + ncp;
    }
    
    if (ncp < 0 || a <= 0 || b <= 0) {
        return InfNaN::nan();
    }
    
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(ncp)) {
        return InfNaN::nan();
    }

    if (x < 0 || x > 1) {
        return give_log ? InfNaN::neginf() : 0.0;
    }
    
    if (ncp == 0) {
        return Beta::pdf(x, a, b, give_log);
    }

    // New algorithm, starting with *largest* term :
    int kMax = 0;
    double ncp2 = 0.5 * ncp;
    double dx2 = ncp2 * x;
    double d = (dx2 - a - 1) / 2;
    double D = d * d + dx2 * (a + b) - a;
    if (D <= 0) {
        kMax = 0;
    }
    else {
        D = ceil(d + sqrt(D));
        kMax = (D > 0) ? (int)D : 0;
    }

    long double term = Beta::pdf(x, a + kMax, b, true);
    long double p_k = Base::poissonPdfRaw(kMax, ncp2, true);
	if (x == 0.0 || !std::isfinite(term) || !std::isfinite((double)p_k)) {
		double tmp = (double)(p_k + term);
		return give_log ? tmp : exp(tmp);
	}

    p_k += term;
	long double sum = term = 1.0;
	long double q = 0.0;
    double k = kMax;
    while (k > 0 && term > sum * eps) {
        k--;
        q = (k + 1) * (k + a) / (k + a + b) / dx2;
        term *= q;
        sum += term;
    }
    
    term = 1.;
    k = kMax;
    do {
        q = dx2 * (k + a + b) / (k + a) / (k + 1);
        k++;
        term *= q;
        sum += term;
    } while (term > sum * eps);

	double tmp = (double)(p_k + logl(sum));
	return give_log ? tmp : exp(tmp);
}

double NonCentralBeta::cdf(double x, double a, double b, double ncp,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(a) || std::isnan(b) || std::isnan(ncp))
		return x + a + b + ncp;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (x <= 0.0) return R_DT_0;
	if (x >= 1.0) return R_DT_1;
	return Base::nonCentralBetaCdf2(x, 1 - x, a, b, ncp, lower_tail, log_p);
}

double NonCentralBeta::quantile(double p, double a, double b, double ncp,
	bool lower_tail, bool log_p)
{
	constexpr double accu = 1e-15;
	constexpr double Eps = 1e-14; /* must be > accu */

	if (std::isnan(p) || std::isnan(a) || std::isnan(b) || std::isnan(ncp))
		return p + a + b + ncp;
	if (!std::isfinite(a))
		return InfNaN::nan();

	if (ncp < 0. || a <= 0. || b <= 0.)
		return InfNaN::nan();

	//R_Q_P01_boundaries(p, 0, 1);
	if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? 1.0 : 0.0;
		if (p == InfNaN::neginf())
			return lower_tail ? 0.0 : 1.0;
	}
	else {
		if (p < 0 || p > 1)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? 0.0 : 1.0;
		if (p == 1)
			return lower_tail ? 1.0 : 0.0;
	}

	const double tmp = lower_tail ? (p) : (0.5 - p +0.5);
	p = log_p ? (lower_tail ? exp(p) : -expm1(p)) : tmp;

	if (p > 1 - DBL_EPSILON)
		return 1.0;

	double ux = 0.0;
	double lx = 0.0;
	double nx = 0.0;
	double pp = std::min(1.0 - DBL_EPSILON, p * (1 + Eps));
	for (ux = 0.5;
		ux < 1 - DBL_EPSILON && cdf(ux, a, b, ncp, true, false) < pp;
		ux = 0.5 * (1 + ux));

	pp = p * (1 - Eps);
	for (lx = 0.5;
		lx > DBL_MIN && cdf(lx, a, b, ncp, true, false) > pp;
		lx *= 0.5);

	do {
		nx = 0.5 * (lx + ux);
		if (cdf(nx, a, b, ncp, true, false) > p)
			ux = nx;
		else
			lx = nx;
	} while ((ux - lx) / nx > accu);

	return 0.5 * (ux + lx);
}

double NonCentralBeta::rand(const double shape1, const double shape2,
	const double ncp)
{
    double x = NonCentralChisq::rand(2.0 * shape1, ncp);
    return x / (x + Chisq::rand(2.0 * shape2));
}