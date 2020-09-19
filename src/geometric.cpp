/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 1999-2016  The R Core Team
 * Copyright (C) 2004-2016 The R Foundation
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

double Geom::pdf(double x, double p, bool give_log)
{
    if (std::isnan(x) || std::isnan(p)) return x + p;

    if (p <= 0 || p > 1) return InfNaN::nan();

    bool isNonInt = fabs((x)- nearbyint(x)) > 1e-7 * std::max(1.0, fabs(x));
    if(isNonInt) {
        std::cout << "Warning: non-integer x = " << x << std::endl;
	    return give_log ? InfNaN::neginf() : 0.0;
    }

    if (x < 0 || !std::isfinite(x) || p == 0)
        return give_log ? InfNaN::neginf() : 0.0;
    x = nearbyint(x);

    /* prob = (1-p)^x, stable for small p */
    double prob = Base::binomialPdfRaw(0.0, x, p, 1 - p, give_log);

    return (give_log) ? log(p) + prob : p * prob;
}

double Geom::cdf(double x, double p, bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(p))
		return x + p;

	if (p <= 0 || p > 1) return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (x < 0.) return R_DT_0;
	if (!std::isfinite(x)) return R_DT_1;
	x = floor(x + 1e-7);

	if (p == 1.0) { /* we cannot assume IEEE */
		x = lower_tail ? 1 : 0;
		return log_p ? log(x) : x;
	}

	x = log1p(-p) * (x + 1);
	double t1 = log_p ? x : log(x);
	double t2 = x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x));
	double t3 = log_p ? t2 : log1p(-x);
	if (log_p)
		return lower_tail ? t3 : t1;
	else
		return lower_tail ? -expm1(x) : exp(x);
}

double Geom::quantile(double p, double prob, bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(prob))
		return p + prob;

	if (prob <= 0 || prob > 1) return InfNaN::nan();

	if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
		return InfNaN::nan();
	if (prob == 1)
		return 0;

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

	// add a fuzz to ensure left continuity, but value must be >= 0
	double t1 = log_p ? p : log(p);
	double t2 = p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p));
	double t3 = log_p ? t2 : log1p(-p);
	double t4 = lower_tail ? t3 : t1;
	double t5 = ceil(t4 / log1p(-prob) - 1.0 - 1e-12);
	return std::max(0.0, t5);
}

int Geom::rand(const double p)
{
    if (!std::isfinite(p) || p <= 0.0 || p > 1.0) {
    	return InfNaN::nan();
    }

    std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::geometric_distribution<int> u(p);
	return u(e);
}