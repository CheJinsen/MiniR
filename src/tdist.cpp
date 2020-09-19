/*
 * This file is part of MiniR.
 * AUTHOR
 *  Catherine Loader, catherine@research.bell-labs.com.
 *  October 23, 2000.
 *
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 2000-2015 The R Core Team
 * Copyright (C) 2003-2013 The R Foundation
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

double Tdist::pdf(double x, double n, bool give_log)
{
	if (std::isnan(x) || std::isnan(n)) {
		return x + n;
	}

	if (n <= 0) {
		return InfNaN::nan();
	}
	
	if (!std::isfinite(x)) {
		return give_log ? InfNaN::neginf() : 0.0;
	}
	
	if (!std::isfinite(n)) {
		return Normal::pdf(x, 0.0, 1.0, give_log);
	}

	double u = 0.0;
	double t = -Base::bd0(n / 2.0, (n + 1) / 2.0) + Base::stirlerr((n + 1) / 2.0) - Base::stirlerr(n / 2.0);
	double x2n = x * x / n;
	double ax = 0.0;
	double l_x2n = 0.0;
	bool lrg_x2n = (x2n > 1.0 / DBL_EPSILON);

	if (lrg_x2n) {
		ax = fabs(x);
		l_x2n = log(ax) - log(n) / 2.0;
		u = n * l_x2n;
	}
	else if (x2n > 0.2) {
		l_x2n = log(1 + x2n) / 2.0;
		u = n * l_x2n;
	}
	else {
		l_x2n = log1p(x2n) / 2.0;
		u = -Base::bd0(n / 2.0, (n + x * x) / 2.0) + x * x / 2.0;
	}

	if (give_log) {
		return t - u - (log(sqrt(2 * M_PI)) + l_x2n);
	}

	double I_sqrt_ = (lrg_x2n ? sqrt(n) / ax : exp(-l_x2n));
	return exp(t - u) * 1.0 / sqrt(2 * M_PI) * I_sqrt_;
}

double Tdist::cdf(double x, double n, bool lower_tail, bool log_p)
{
	double val = 0.0, nx = 0.0;

	if (std::isnan(x) || std::isnan(n)) {
		return x + n;
	}

	if (n <= 0.0) {
		return InfNaN::nan();
	}

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (!std::isfinite(x)) {
		return (x < 0) ? R_DT_0 : R_DT_1;
	}
	
	if (!std::isfinite(n)) {
		return Normal::cdf(x, 0.0, 1.0, lower_tail, log_p);
	}

	nx = 1 + (x / n) * x;
	if (nx > 1e100) {
		double lval;
		lval = -0.5 * n * (2 * log(fabs(x)) - log(n))
			- SpecialFunctions::Beta::lbeta(0.5 * n, 0.5) - log(0.5 * n);
		val = log_p ? lval : exp(lval);
	}
	else {
		val = (n > x * x)
			? Beta::cdf(x * x / (n + x * x), 0.5, n / 2.0, false, log_p)
			: Beta::cdf(1.0 / nx, n / 2.0, 0.5, true, log_p);
	}

	if (x <= 0.0) {
		lower_tail = !lower_tail;
	}

	if (log_p) {
		if (lower_tail) return log1p(-0.5 * exp(val));
		else 			return val - M_LN2;
	}
	else {
		val /= 2.;
		return lower_tail ? (0.5 - val + 0.5) : val;
	}
}

double Tdist::quantile(double p, double ndf, bool lower_tail, bool log_p)
{
	const double eps = 1.e-12;

	double P = 0.0, q = 0.0;
	if (std::isnan(p) || std::isnan(ndf)) {
		return p + ndf;
	}

	//R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
	if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? InfNaN::posinf() : InfNaN::neginf();
		if (p == InfNaN::neginf())
			return lower_tail ? InfNaN::neginf() : InfNaN::posinf();
	}
	else {
		if (p < 0 || p > 1)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? InfNaN::neginf() : InfNaN::posinf();
		if (p == 1)
			return lower_tail ? InfNaN::posinf() : InfNaN::neginf();
	}

	if (ndf <= 0) {
		return InfNaN::nan();
	}

	if (ndf < 1) {
		const double accu = 1e-13;
		const double Eps = 1e-11;
		double ux = 0.0, lx = 0.0, nx = 0.0, pp = 0.0;
		int iter = 0;

		p = log_p ? (lower_tail ? exp(p) : - expm1(p))
			: (lower_tail ? (p) : (0.5 - (p)+0.5));

		if (p > 1 - DBL_EPSILON) return InfNaN::posinf();
		pp = std::min(1 - DBL_EPSILON, p * (1 + Eps));
		for (ux = 1.; ux < DBL_MAX && cdf(ux, ndf, true, false) < pp; ux *= 2);
		pp = p * (1 - Eps);
		for (lx = -1.; lx > -DBL_MAX && cdf(lx, ndf, true, false) > pp; lx *= 2);

		do {
			nx = 0.5 * (lx + ux);
			if (cdf(nx, ndf, true, false) > p) ux = nx; else lx = nx;
		} while ((ux - lx) / fabs(nx) > accu && ++iter < 1000);

		if (iter >= 1000)
			std::cout << "Full precision may not have been "
				<< "achieved in Tdist quantile()" << std::endl;

		return 0.5 * (lx + ux);
	}

	if (ndf > 1e20) {
		return Normal::quantile(p, 0.0, 1.0, lower_tail, log_p);
	}

	P = log_p ? exp(p) : p;
	bool neg = (!lower_tail || P < 0.5) && (lower_tail || P > 0.5);
	bool is_neg_lower = (lower_tail == neg);
	
	if (neg) {
		P = 2 * (log_p ? (lower_tail ? P : -expm1(p))
			: (lower_tail ? p : (0.5 - p + 0.5)));
	}
	else {
		P = 2 * (log_p ? (lower_tail ? -expm1(p) : P)
			: (lower_tail ? (0.5 - p + 0.5) : p));
	}

	if (fabs(ndf - 2) < eps) {
		if (P > DBL_MIN) {
			if (3 * P < DBL_EPSILON)
				q = 1 / sqrt(P);
			else if (P > 0.9)
				q = (1 - P) * sqrt(2 / (P * (2 - P)));
			else
				q = sqrt(2 / (P * (2 - P)) - 2);
		}
		else {
			if (log_p)
				q = is_neg_lower ? exp(-p / 2) / sqrt(2.0) : 1 / sqrt(-expm1(p));
			else
				q = InfNaN::posinf();
		}
	}
	else if (ndf < 1 + eps) {
		if (P == 1.0)
			q = 0;
		else if (P > 0)
			q = 1 / Base::tanpi(P / 2.);

		else {
			if (log_p)
				q = is_neg_lower ? M_1_PI * exp(-p) : -1. / (M_PI * expm1(p));
			else
				q = InfNaN::posinf();
		}
	}
	else {		/*-- usual case;  including, e.g.,  df = 1.1 */
		double x = 0., y, log_P2 = 0./* -Wall */,
			a = 1 / (ndf - 0.5),
			b = 48 / (a * a),
			c = ((20700 * a / b - 98) * a - 16) * a + 96.36,
			d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;

		bool P_ok1 = P > DBL_MIN || !log_p, P_ok = P_ok1;
		if (P_ok1) {
			y = pow(d * P, 2.0 / ndf);
			P_ok = (y >= DBL_EPSILON);
		}

		double temp1 = log_p ? p : log(p);
		double temp2 = p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p));
		double temp3 = log_p ? temp2 : log1p(-p);

		if (!P_ok) {
			log_P2 = is_neg_lower ? temp1 : temp3;
			x = (log(d) + M_LN2 + log_P2) / ndf;
			y = exp(2 * x);
		}

		if ((ndf < 2.1 && P > 0.5) || y > 0.05 + a) {
			if (P_ok)
				x = Normal::quantile(0.5 * P, 0.0, 1.0, true, false);
			else
				x = Normal::quantile(log_P2, 0., 1., lower_tail, true);

			y = x * x;
			if (ndf < 5)
				c += 0.3 * (ndf - 4.5) * (x + 0.6);
			c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
			y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c
				- y - 3) / b + 1) * x;
			y = expm1(a * y * y);
			q = sqrt(ndf * y);
		}
		else if (!P_ok && x < -M_LN2 * DBL_MANT_DIG) {
			q = sqrt(ndf) * exp(-x);
		}
		else {
			y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
				* (ndf + 2) * 3) + 0.5 / (ndf + 4))
				* y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
			q = sqrt(ndf * y);
		}

		if (P_ok1) {
			int it = 0;
			while (it++ < 10 && (y = pdf(q, ndf, false)) > 0 &&
				std::isfinite(x = (cdf(q, ndf, false, false) - P / 2) / y) &&
				fabs(x) > 1e-14 * fabs(q))
					q += x * (1. + x * q * (ndf + 1) / (2 * (q * q + ndf)));
		}
	}
	if (neg) {
		q = -q;
	}
	
	return q;
}

double Tdist::rand(const double df)
{
	if (std::isnan(df) || df <= 0.0) {
		return InfNaN::nan();
	}

	std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::student_t_distribution<double> u(df);
	return u(e);
}