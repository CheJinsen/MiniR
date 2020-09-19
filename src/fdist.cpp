/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 1999-2015 The R Core Team
 * Copyright (C) 2004-2005 The R Foundation
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

double Fdist::pdf(double x, double m, double n, bool give_log)
{	
	double dens = 0.0;
	if (std::isnan(x) || std::isnan(m) || std::isnan(n))
		return x + m + n;
	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
	double R_D__1 = give_log ? 0.0 : 1.0;

	if (m <= 0 || n <= 0) return InfNaN::nan();
	if (x < 0.)  return R_D__0;
	if (x == 0.) return m > 2 ? R_D__0 : (m == 2 ? R_D__1 : InfNaN::posinf());
	if (!std::isfinite(m) && !std::isfinite(n)) {
		if (x == 1.0) return InfNaN::posinf(); else return R_D__0;
	}
	if (!std::isfinite(n))
		return(Gamma::pdf(x, m / 2, 2.0 / m, give_log));
	if (m > 1e14) {
		dens = Gamma::pdf(1.0 / x, n / 2, 2.0 / n, give_log);
		return give_log ? dens - 2 * log(x) : dens / (x * x);
	}

	double f = 1.0 / (n + x * m);
	double q = n * f;
	double p = x * m * f;

	if (m >= 2) {
		f = m * q / 2;
		dens = Base::binomialPdfRaw((m - 2) / 2, (m + n - 2) / 2, p, q, give_log);
	}
	else {
		f = m * m * q / (2 * p * (m + n));
		dens = Base::binomialPdfRaw(m / 2, (m + n) / 2, p, q, give_log);
	}
	return(give_log ? log(f) + dens : f * dens);
}

double Fdist::cdf(double x, double df1, double df2,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(df1) || std::isnan(df2))
		return x + df2 + df1;
	if (df1 <= 0.0 || df2 <= 0.0) return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if(x <= 0.0) return R_DT_0;
	if (x >= InfNaN::posinf()) return R_DT_1;

	if (df2 == InfNaN::posinf()) {
		if (df1 == InfNaN::posinf()) {
			if (x < 1.0) return R_DT_0;
			if (x == 1.0) return (log_p ? -M_LN2 : 0.5);
			if (x > 1.0) return R_DT_1;
		}
		return Chisq::cdf(x * df1, df1, lower_tail, log_p);
	}

	if (df1 == InfNaN::posinf())
		return Chisq::cdf(df2 / x, df2, !lower_tail, log_p);

	if (df1 * x > df2)
		x = Beta::cdf(df2 / (df2 + df1 * x), df2 / 2., df1 / 2.,
			!lower_tail, log_p);
	else
		x = Beta::cdf(df1 * x / (df2 + df1 * x), df1 / 2., df2 / 2.,
			lower_tail, log_p);

	return !std::isnan(x) ? x : InfNaN::nan();
}

double Fdist::quantile(double p, double df1, double df2,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(df1) || std::isnan(df2))
		return p + df1 + df2;

	if (df1 <= 0. || df2 <= 0.) return InfNaN::nan();

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

	if (df1 <= df2 && df2 > 4e5) {
		if (!std::isfinite(df1))
			return 1.0;
		return Chisq::quantile(p, df1, lower_tail, log_p) / df1;
	}
	if (df1 > 4e5) {
		return df2 / Chisq::quantile(p, df2, !lower_tail, log_p);
	}

	p = (1.0 / Beta::quantile(p, df2 / 2, df1 / 2, !lower_tail, log_p) - 1.0) * (df2 / df1);
	return !std::isnan(p) ? p : InfNaN::nan();
}

double Fdist::rand(const double n1, const double n2)
{
    if (std::isnan(n1) || std::isnan(n2) || n1 <= 0.0 || n2 <= 0.0) {
		return InfNaN::nan();
    }

    std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::fisher_f_distribution<double> u(n1, n2);
	return u(e);
}