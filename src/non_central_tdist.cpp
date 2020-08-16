/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Claus Ekstr鴐, ekstrom@dina.kvl.dk
 *    July 15, 2003.
 *
 * Copyright (C) 2003-2015 The R Foundation
 * Copyright (C) 1998-2015 The R Core Team
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

double NonCentralTdist::pdf(double x, double df, double ncp, bool give_log)
{
    if (std::isnan(x) || std::isnan(df)) {
		return x + df;
    }

    /* If non-positive df then error */
    if (df <= 0.0) return InfNaN::nan();
    if (ncp == 0.0) return Tdist::pdf(x, df, give_log);

    if (!std::isfinite(x))
		return give_log ? InfNaN::neginf() : 0.0;;

    if (!std::isfinite(df) || df > 1e8)
		return Normal::pdf(x, ncp, 1.0, give_log);

	double u = 0.0;
    if (fabs(x) > sqrt(df * DBL_EPSILON)) {
		u = log(df) - log(fabs(x)) +
		    log(fabs(cdf(x * sqrt((df + 2) / df), df + 2, ncp, true, false) -
			cdf(x, df, ncp, true, false)));
    }
    else {
		u = SpecialFunctions::Gamma::lgammafn((df + 1) / 2) -
			SpecialFunctions::Gamma::lgammafn(df / 2) -
			(log(sqrt(M_PI)) + 0.5 * (log(df) + ncp * ncp));
    }

    return give_log ? u : exp(u);
}

double NonCentralTdist::cdf(double t, double df, double ncp, bool lower_tail, bool log_p)
{
    double albeta = 0.0, a = 0.0, b = 0.0, del = 0.0;
    double errbd = 0.0, lambda = 0.0, tt = 0.0;
    long double geven = 0.0, godd = 0.0, p = 0.0, q = 0.0;
    long double s = 0.0, tnc = 0.0, xeven = 0.0, xodd = 0.0;
    bool negdel = false;

    constexpr int itrmax = 1000;
    constexpr double errmax = 1.e-12;

    if (df <= 0.0) 	return InfNaN::nan();
    if (ncp == 0.0) return Tdist::cdf(t, df, lower_tail, log_p);

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (!std::isfinite(t)) {
		return (t < 0) ? R_DT_0 : R_DT_1;
    }
    if (t >= 0.0) {
		negdel = false; tt = t;	 del = ncp;
    }
    else {
		if (ncp > 40 && (!log_p || !lower_tail))
			return R_DT_0;
		negdel = true;	tt = -t; del = -ncp;
    }

    if (df > 4e5 || del * del > 2 * M_LN2 * (-(DBL_MIN_EXP))) {
		s = 1.0 / (4.0 * df);
		return Normal::cdf((double)(tt * (1.0 - s)), del,
			     sqrt((double)(1.0 + tt * tt * 2.0 * s)),
			     lower_tail != negdel, log_p);
    }

    double x = t * t;
    double rxb = df/(x + df);
    x = x / (x + df);
    if (x > 0.0) {
		lambda = del * del;
		p = 0.5 * exp(-0.5 * lambda);

		if (p == 0.0) {
		    std::cout << "Warning: Underflow occurred in cdf()." << std::endl;
		    std::cout << "Warning: Value out of range in cdf()." << std::endl;
		    return R_DT_0;
		}
	
		q = sqrt(2.0 / M_PI) * p * del;
		s = 0.5 - p;
		
		if (s < 1e-7)
		    s = -0.5 * expm1(-0.5 * lambda);
		a = 0.5;
		b = 0.5 * df;

		rxb = pow(rxb, b);
		albeta = log(sqrt(M_PI)) + SpecialFunctions::Gamma::lgammafn(b) -
			SpecialFunctions::Gamma::lgammafn(0.5 + b);
		xodd = Beta::cdf(x, a, b, true, false);
		godd = 2.0 * rxb * exp(a * log(x) - albeta);
		tnc = b * x;
		xeven = (tnc < DBL_EPSILON) ? tnc : 1.0 - rxb;
		geven = tnc * rxb;
		tnc = p * xodd + q * xeven;

		 // repeat until convergence or iteration limit 
		for (int it = 1; it <= itrmax; it++) {
		    a += 1.0;
		    xodd  -= godd;
		    xeven -= geven;
		    godd  *= x * (a + b - 1.0) / a;
		    geven *= x * (a + b - 0.5) / (a + 0.5);
		    p *= lambda / (2 * it);
		    q *= lambda / (2 * it + 1);
		    tnc += p * xodd + q * xeven;
		    s -= p;

		    if (s < -1.e-10) {
				std::cout << "Warning: Full precision may not have "
					<< "been achieved in cdf()" << std::endl;
				goto finis;
		    }

		    if (s <= 0 && it > 1)
		    	goto finis;

		    errbd = 2.0 * s * (xodd - godd);
		    if (fabs(errbd) < errmax)
		    	goto finis;
		}
		std::cout << "Warning: convergence failed in cdf()." << std::endl;
    }
    else { /* x = t = 0 */
		tnc = 0.0;
    }

 finis:
    tnc += Normal::cdf(-del, 0.0, 1.0, true, false);

    lower_tail = lower_tail != negdel;
    if (tnc > 1 - 1e-10 && lower_tail) {
		std::cout << "Warning: Full precision may not have "
					<< "been achieved in cdf()" << std::endl;
    }

    double value = std::min((double)tnc, 1.0);
    double temp1 = log_p ? log(value) : value;
    double temp2 = log_p ? log1p(-value) : (0.5 - value + 0.5);
    return lower_tail ? temp1 : temp2;
}

double NonCentralTdist::quantile(double p, double df, double ncp, bool lower_tail, bool log_p)
{
    constexpr double accu = 1e-13;
    constexpr double Eps = 1e-11; /* must be > accu */

    double ux = 0.0, lx = 0.0, nx = 0.0;

    if (std::isnan(p) || std::isnan(df) || std::isnan(ncp))
		return p + df + ncp;

    if (df <= 0.0) return InfNaN::nan();

    if(ncp == 0.0 && df >= 1.0)
    	return Tdist::quantile(p, df, lower_tail, log_p);

    // R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
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

    if (!std::isfinite(df)) // df = Inf ==> limit N(ncp,1)
		return Normal::quantile(p, ncp, 1.0, lower_tail, log_p);

	double temp = lower_tail ? p : (0.5 - p + 0.5);
    p = log_p ? (lower_tail ? exp(p) : - expm1(p)) : temp;


    if(p > 1 - DBL_EPSILON)
    	return InfNaN::posinf();

    double pp = std::min(1.0 - DBL_EPSILON, p * (1.0 + Eps));

    for(ux = std::max(1.0, ncp);
    	ux < DBL_MAX && cdf(ux, df, ncp, true, false) < pp; ux *= 2)
    	;

    pp = p * (1 - Eps);
    for(lx = std::min(-1.0, -ncp);
		lx > -DBL_MAX && cdf(lx, df, ncp, true, false) > pp; lx *= 2)
		;

    do {
		nx = 0.5 * (lx + ux); // could be zero
		if (cdf(nx, df, ncp, true, false) > p)
			ux = nx;
		else
			lx = nx;
    } while ((ux - lx) > accu * std::max(fabs(lx), fabs(ux)));

    return 0.5 * (lx + ux);
}

double NonCentralTdist::rand(const double df, const double ncp)
{
	return Normal::rand(ncp) / sqrt(Chisq::rand(df) / df);
}