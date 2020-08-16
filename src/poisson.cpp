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

double Poisson::pdfRaw(double x, double lambda, bool give_log)
{
    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
    double R_D__1 = give_log ? 0.0 : 1.0;

    if (lambda == 0) return((x == 0) ? R_D__1 : R_D__0);
    if (!std::isfinite(lambda)) return R_D__0;
    if (x < 0) return(R_D__0);
	if (x <= lambda * DBL_MIN) return give_log ? -lambda : exp(-lambda);
    if (lambda < x * DBL_MIN) {
        if (!std::isfinite(x))
            return R_D__0;
        double temp = -lambda + x * log(lambda) - SpecialFunctions::Gamma::lgammafn(x + 1);
		return give_log ? temp : exp(temp);
    }
	double f = 2.0 * M_PI * x;
	x = -stirlerr(x) - bd0(x, lambda);
	return give_log ? -0.5 * log(f) + (x) : exp(x) / sqrt(f);
}

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
    return pdfRaw(x, lambda, give_log);
}

double Poisson::bd0(const double x, const double np)
{
	double ej = 0.0, s = 0.0, s1 = 0.0, v = 0.0;

	if (!std::isfinite(x) || !std::isfinite(np) || np == 0.0)
		return InfNaN::nan();

	if (fabs(x - np) < 0.1 * (x + np)) {
		v = (x - np) / (x + np);  // might underflow to 0
		s = (x - np) * v;
		if (fabs(s) < DBL_MIN) return s;
		ej = 2 * x * v;
		v = v * v;
		for (int j = 1; j < 1000; j++) {
			ej *= v;
			s1 = s + ej / ((j * 2.0) + 1);
			if (s1 == s)
				return s1;
			s = s1;
		}
	}
	return(x * log(x / np) + np - x);
}

double Poisson::stirlerr(const double n)
{
	constexpr auto S0 = 0.083333333333333333333;       /* 1/12 */;
	constexpr auto S1 = 0.00277777777777777777778;     /* 1/360 */;
	constexpr auto S2 = 0.00079365079365079365079365;  /* 1/1260 */;
	constexpr auto S3 = 0.000595238095238095238095238; /* 1/1680 */;
	constexpr auto S4 = 0.0008417508417508417508417508;/* 1/1188 */;

	/*
	  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
	*/
	const double sferr_halves[31] = {
		0.0, /* n=0 - wrong, place holder only */
		0.1534264097200273452913848,  /* 0.5 */
		0.0810614667953272582196702,  /* 1.0 */
		0.0548141210519176538961390,  /* 1.5 */
		0.0413406959554092940938221,  /* 2.0 */
		0.03316287351993628748511048, /* 2.5 */
		0.02767792568499833914878929, /* 3.0 */
		0.02374616365629749597132920, /* 3.5 */
		0.02079067210376509311152277, /* 4.0 */
		0.01848845053267318523077934, /* 4.5 */
		0.01664469118982119216319487, /* 5.0 */
		0.01513497322191737887351255, /* 5.5 */
		0.01387612882307074799874573, /* 6.0 */
		0.01281046524292022692424986, /* 6.5 */
		0.01189670994589177009505572, /* 7.0 */
		0.01110455975820691732662991, /* 7.5 */
		0.010411265261972096497478567, /* 8.0 */
		0.009799416126158803298389475, /* 8.5 */
		0.009255462182712732917728637, /* 9.0 */
		0.008768700134139385462952823, /* 9.5 */
		0.008330563433362871256469318, /* 10.0 */
		0.007934114564314020547248100, /* 10.5 */
		0.007573675487951840794972024, /* 11.0 */
		0.007244554301320383179543912, /* 11.5 */
		0.006942840107209529865664152, /* 12.0 */
		0.006665247032707682442354394, /* 12.5 */
		0.006408994188004207068439631, /* 13.0 */
		0.006171712263039457647532867, /* 13.5 */
		0.005951370112758847735624416, /* 14.0 */
		0.005746216513010115682023589, /* 14.5 */
		0.005554733551962801371038690  /* 15.0 */
	};
	double nn = 0.0;

	if (n <= 15.0) {
		nn = n + n;
		if (nn == (int)nn) return(sferr_halves[(int)nn]);
		return SpecialFunctions::Gamma::lgammafn(n + 1.0) - (n + 0.5) * log(n) + n - log(sqrt(2 * M_PI));
	}

	nn = n * n;
	if (n > 500) return (S0 - S1 / nn) / n;
	if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
	if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
	/* 15 < n <= 35 : */
	return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
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