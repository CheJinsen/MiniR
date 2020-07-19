/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 1999-2019 The R Core Team
 * Copyright (C) 2004-2019 The R Foundation
 * Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
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

double Gamma::pdf(double x, double shape, double scale, bool give_log)
{
	if (std::isnan(x) || std::isnan(shape) || std::isnan(scale))
		return x + shape + scale;

	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
	if (shape < 0 || scale <= 0) InfNaN::nan();
	if (x < 0)
		return R_D__0;
	if (shape == 0) /* point mass at 0 */
		return (x == 0) ? InfNaN::posinf() : R_D__0;
	if (x == 0) {
		if (shape < 1) return InfNaN::posinf();
		if (shape > 1) return R_D__0;
		return give_log ? -log(scale) : 1 / scale;
	}

	double pr = 0.0;
	if (shape < 1) {
		pr = poissonPdfRaw(shape, x / scale, give_log);
		return (
			give_log/* NB: currently *always*  shape/x > 0  if shape < 1:
				 * -- overflow to Inf happens, but underflow to 0 does NOT : */
			? pr + (std::isfinite(shape / x)
				? log(shape / x)
				: log(shape) - log(x))
			: pr * shape / x);
	}

	pr = poissonPdfRaw(shape - 1, x / scale, give_log);
	return give_log ? pr - log(scale) : pr / scale;
}

double Gamma::poissonPdfRaw(double x, double lambda, bool give_log)
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

    double f = M_2PI * x;
    x = -stirlerr(x) - bd0(x, lambda);
    return give_log ? -0.5 * log(f) + (x) : exp(x) / sqrt(f);
}

double Gamma::bd0(const double x, const double np)
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

double Gamma::stirlerr(const double n)
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
		return(SpecialFunctions::Gamma::lgammafn(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI);
	}

	nn = n * n;
	if (n > 500) return (S0 - S1 / nn) / n;
	if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
	if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
	/* 15 < n <= 35 : */
	return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

double Gamma::logcf(double x, double i, double d, double eps)
{
	double c1 = 2 * d;
	double c2 = i + d;
	double c4 = c2 + d;
	double a1 = c2;
	double b1 = i * (c2 - i * x);
	double b2 = d * d * x;
	double a2 = c4 * c2 - b2;

	b2 = c4 * b1 - i * b2;
	double scalefactor = sqr(sqr(sqr(4294967296.0)));

	while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
		double c3 = c2 * c2 * x;
		c2 += d;
		c4 += d;
		a1 = c4 * a2 - c3 * a1;
		b1 = c4 * b2 - c3 * b1;

		c3 = c1 * c1 * x;
		c1 += d;
		c4 += d;
		a2 = c4 * a1 - c3 * a2;
		b2 = c4 * b1 - c3 * b2;

		if (fabs(b2) > scalefactor) {
			a1 /= scalefactor;
			b1 /= scalefactor;
			a2 /= scalefactor;
			b2 /= scalefactor;
		}
		else if (fabs(b2) < 1 / scalefactor) {
			a1 *= scalefactor;
			b1 *= scalefactor;
			a2 *= scalefactor;
			b2 *= scalefactor;
		}
	}
	return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double Gamma::log1pmx(double x)
{
	const double minLog1Value = -0.79149064;
	if (x > 1 || x < minLog1Value) {
		return log1p(x) - x;
	}
	else {
		double r = x / (2 + x);
		double y = r * r;
		if (fabs(x) < 1e-2) {
			static const double two = 2;
			return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
				two / 3) * y - x);
		}
		else {
			const double tol_logcf = 1e-14;
			return r * (2 * y * logcf(y, 3, 2, tol_logcf) - x);
		}
	}
}

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double Gamma::lgamma1p(double a)
{
	const double eulers_const = 0.5772156649015328606065120900824024;

	/* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
	const int N = 40;
	const double coeffs[40] = {
		0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
		0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
		0.2058080842778454787900092413529198e-1,
		0.7385551028673985266273097291406834e-2,
		0.2890510330741523285752988298486755e-2,
		0.1192753911703260977113935692828109e-2,
		0.5096695247430424223356548135815582e-3,
		0.2231547584535793797614188036013401e-3,
		0.9945751278180853371459589003190170e-4,
		0.4492623673813314170020750240635786e-4,
		0.2050721277567069155316650397830591e-4,
		0.9439488275268395903987425104415055e-5,
		0.4374866789907487804181793223952411e-5,
		0.2039215753801366236781900709670839e-5,
		0.9551412130407419832857179772951265e-6,
		0.4492469198764566043294290331193655e-6,
		0.2120718480555466586923135901077628e-6,
		0.1004322482396809960872083050053344e-6,
		0.4769810169363980565760193417246730e-7,
		0.2271109460894316491031998116062124e-7,
		0.1083865921489695409107491757968159e-7,
		0.5183475041970046655121248647057669e-8,
		0.2483674543802478317185008663991718e-8,
		0.1192140140586091207442548202774640e-8,
		0.5731367241678862013330194857961011e-9,
		0.2759522885124233145178149692816341e-9,
		0.1330476437424448948149715720858008e-9,
		0.6422964563838100022082448087644648e-10,
		0.3104424774732227276239215783404066e-10,
		0.1502138408075414217093301048780668e-10,
		0.7275974480239079662504549924814047e-11,
		0.3527742476575915083615072228655483e-11,
		0.1711991790559617908601084114443031e-11,
		0.8315385841420284819798357793954418e-12,
		0.4042200525289440065536008957032895e-12,
		0.1966475631096616490411045679010286e-12,
		0.9573630387838555763782200936508615e-13,
		0.4664076026428374224576492565974577e-13,
		0.2273736960065972320633279596737272e-13,
		0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
	};

	const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
	const double tol_logcf = 1e-14;

	if (fabs(a) >= 0.5)
		return SpecialFunctions::Gamma::lgammafn(a + 1);

	double lgam = c * logcf(-a / 2, 1.0 * N + 2, 1, tol_logcf);
	for (int i = N - 1; i >= 0; i--)
		lgam = coeffs[i] - a * lgam;

	return (a * lgam - eulers_const) * a - log1pmx(a);
} /* lgamma1p */

double Gamma::logspaceAdd(double logx, double logy)
{
	return std::max(logx, logy) + log1p(exp(-fabs(logx - logy)));
}

double Gamma::logspaceSub(double logx, double logy)
{
	double x = logy - logx;
	double temp = x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x));
	return logx + temp;
}

double Gamma::logspaceSum(const double* logx, int n)
{
	if (n == 0) return InfNaN::neginf(); // = log( sum(<empty>) )
	if (n == 1) return logx[0];
	if (n == 2) return logspaceAdd(logx[0], logx[1]);
	
	double Mx = logx[0];
	for (int i = 1; i < n; i++)
		if (Mx < logx[i])
			Mx = logx[i];
	long double s = (long double)0.0;
	for (int i = 0; i < n; i++)
		s += expl(logx[i] - Mx);
	return Mx + (double)logl(s);
}

double Gamma::poissonPdfWrap(double x_plus_1, double lambda, bool give_log)
{
	const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;	// 3.196577e18
	
	if (!std::isfinite(lambda))
		return give_log ? InfNaN::neginf() : 0.0;
	if (x_plus_1 > 1)
		return poissonPdfRaw(x_plus_1 - 1, lambda, give_log);

	if (lambda > fabs(x_plus_1 - 1) * M_cutoff) {
		double x = -lambda - SpecialFunctions::Gamma::lgammafn(x_plus_1);
		return give_log ? x : exp(x);
	}
	else {
		double d = poissonPdfRaw(x_plus_1, lambda, give_log);

		return give_log
			? d + log(x_plus_1 / lambda)
			: d * (x_plus_1 / lambda);
	}
}

double Gamma::cdfSmallx(double x, double alph, bool lower_tail, bool log_p)
{
	double sum = 0.0, c = alph, n = 0.0, term = 0.0;

	do {
		n++;
		c *= -x / n;
		term = c / (alph + n);
		sum += term;
	} while (fabs(term) > DBL_EPSILON * fabs(sum));

	if (lower_tail) {
		double f1 = log_p ? log1p(sum) : 1 + sum;
		double f2 = 0.0;
		if (alph > 1) {
			f2 = poissonPdfRaw(alph, x, log_p);
			f2 = log_p ? f2 + x : f2 * exp(x);
		}
		else if (log_p)
			f2 = alph * log(x) - lgamma1p(alph);
		else
			f2 = pow(x, alph) / exp(lgamma1p(alph));

		return log_p ? f1 + f2 : f1 * f2;
	}
	else {
		double lf2 = alph * log(x) - lgamma1p(alph);

		if (log_p) {
			double temp = log1p(sum) + lf2;
			return temp > -M_LN2 ? log(-expm1(temp)) : log1p(-exp(temp));
		}
		else {
			double f1m1 = sum;
			double f2m1 = expm1(lf2);
			return -(f1m1 + f2m1 + f1m1 * f2m1);
		}
	}
} /* pgamma_smallx() */

double Gamma::pdUpperSeries(double x, double y, bool log_p)
{
	double term = x / y;
	double sum = term;

	do {
		y++;
		term *= x / y;
		sum += term;
	} while (term > sum * DBL_EPSILON);

	return log_p ? log(sum) : sum;
}

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
double Gamma::pdLowercf(double y, double d)
{
	double f = 0.0;
	double scalefactor = sqr(sqr(sqr(4294967296.0)));
	constexpr auto max_it = 200000;

	if (y == 0) return 0;

	double f0 = y / d;
	if (fabs(y - 1) < fabs(d) * DBL_EPSILON) {
		return f0;
	}

	if (f0 > 1.0) f0 = 1.0;
	double c2 = y;
	double c4 = d; /* original (y,d), *not* potentially scaled ones!*/

	double a1 = 0.0, b1 = 1.0;
	double a2 = y, b2 = d;

	while (b2 > scalefactor) {
		a1 /= scalefactor;
		b1 /= scalefactor;
		a2 /= scalefactor;
		b2 /= scalefactor;
	}

	double i = 0.0; 
	double of = -1.0; /* far away */
	while (i < max_it) {

		i++;	c2--;	double c3 = i * c2;	c4 += 2;
		/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
		a1 = c4 * a2 + c3 * a1;
		b1 = c4 * b2 + c3 * b1;

		i++;	c2--;	c3 = i * c2;	c4 += 2;
		/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
		a2 = c4 * a1 + c3 * a2;
		b2 = c4 * b1 + c3 * b2;

		if (b2 > scalefactor) {
			a1 /= scalefactor;
			b1 /= scalefactor;
			a2 /= scalefactor;
			b2 /= scalefactor;
		}

		if (b2 != 0) {
			f = a2 / b2;
			/* convergence check: relative; "absolute" for very small f : */
			if (fabs(f - of) <= DBL_EPSILON * std::max(f0, fabs(f))) {
				return f;
			}
			of = f;
		}
	}

	std::cout << " ** NON-convergence in pgamma()'s pd_lower_cf() f= "
		<< f << std::endl;
	return f;/* should not happen ... */
}

double Gamma::pdLowerSeries(double lambda, double y)
{
	double term = 1.0, sum = 0.0;

	while (y >= 1 && term > sum * DBL_EPSILON) {
		term *= y / lambda;
		sum += term;
		y--;
	}

	if (y != floor(y)) {
		double f = pdLowercf(y, lambda + 1 - y);
		sum += term * f;
	}

	return sum;
}

double Gamma::dpnorm(double x, bool lower_tail, double lp)
{
	if (x < 0) {
		x = -x;
		lower_tail = !lower_tail;
	}

	if (x > 10 && !lower_tail) {
		double term = 1 / x;
		double sum = term;
		double x2 = x * x;
		double i = 1;

		do {
			term *= -i / x2;
			sum += term;
			i += 2;
		} while (fabs(term) > DBL_EPSILON * sum);

		return 1 / sum;
	}
	else {
		double d = Normal::pdf(x, 0.0, 1.0, false);
		return d / exp(lp);
	}
}

double Gamma::poissonCdfAsymp(double x, double lambda, bool lower_tail, bool log_p)
{
	const double coefs_a[8] = {
		-1e99, /* placeholder used for 1-indexing */
		2 / 3.0,
		-4 / 135.0,
		8 / 2835.0,
		16 / 8505.0,
		-8992 / 12629925.0,
		-334144 / 492567075.0,
		698752 / 1477701225.0
	};

	const double coefs_b[8] = {
		-1e99, /* placeholder */
		1 / 12.0,
		1 / 288.0,
		-139 / 51840.0,
		-571 / 2488320.0,
		163879 / 209018880.0,
		5246819 / 75246796800.0,
		-534703531 / 902961561600.0
	};

	double dfm = lambda - x;
	double pt_ = -log1pmx(dfm / x);
	double s2pt = sqrt(2 * x * pt_);
	if (dfm < 0) s2pt = -s2pt;

	double res12 = 0;
	double res1_ig = sqrt(x);
	double res1_term = sqrt(x);
	double res2_ig = s2pt;
	double res2_term = s2pt;
	for (int i = 1; i < 8; i++) {
		res12 += res1_ig * coefs_a[i];
		res12 += res2_ig * coefs_b[i];
		res1_term *= pt_ / i;
		res2_term *= 2.0 * pt_ / (2.0 * i + 1);
		res1_ig = res1_ig / x + res1_term;
		res2_ig = res2_ig / x + res2_term;
	}

	double elfb = x;
	double elfb_term = 1;
	for (int i = 1; i < 8; i++) {
		elfb += elfb_term * coefs_b[i];
		elfb_term /= x;
	}
	if (!lower_tail) elfb = -elfb;

	double f = res12 / elfb;
	double np = Normal::cdf(s2pt, 0.0, 1.0, !lower_tail, log_p);

	if (log_p) {
		double n_d_over_p = dpnorm(s2pt, !lower_tail, np);
		return np + log1p(f * n_d_over_p);
	}
	else {
		double nd = Normal::pdf(s2pt, 0.0, 1.0, log_p);
		return np + f * nd;
	}
}

double Gamma::cdfRaw(double x, double alph, bool lower_tail, bool log_p)
{
	/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

	double res = 0.0;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	
	if (x <= 0.0) return R_DT_0;
	if (x >= InfNaN::posinf()) return R_DT_1;

	if (x < 1) {
		res = cdfSmallx(x, alph, lower_tail, log_p);
	}
	else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
		double sum = pdUpperSeries(x, alph, log_p);
		double d = poissonPdfWrap(alph, x, log_p);

		double temp = d + sum;
		double temp2 = temp > -M_LN2 ? log(-expm1(temp)) : log1p(-exp(temp));
		if (!lower_tail)
			res = log_p ? temp2 : 1 - d * sum;
		else
			res = log_p ? sum + d : sum * d;
	}
	else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
		/* incl. large x compared to alph */
		double sum = 0.0;
		double d = poissonPdfWrap(alph, x, log_p);

		if (alph < 1) {
			if (x * DBL_EPSILON > 1 - alph)
				sum = R_D__1;
			else {
				double f = pdLowercf(alph, x - (alph - 1)) * x / alph;
				/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
				sum = log_p ? log(f) : f;
			}
		}
		else {
			sum = pdLowerSeries(x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
			sum = log_p ? log1p(sum) : 1 + sum;
		}

		if (!lower_tail) {
			res = log_p ? sum + d : sum * d;
		}
		else {
			double temp = d + sum;
			double temp2 = temp > -M_LN2 ? log(-expm1(temp)) : log1p(-exp(temp));
			res = log_p ? temp2 : 1 - d * sum;
		}
	}
	else {
		res = poissonCdfAsymp(alph - 1, x, !lower_tail, log_p);
	}

	if (!log_p && res < DBL_MIN / DBL_EPSILON) {
		return exp(cdfRaw(x, alph, lower_tail, 1));
	}
	else
		return res;
}

double Gamma::cdf(double x, double alph, double scale, bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(alph) || std::isnan(scale))
		return x + alph + scale;

	if (alph < 0. || scale <= 0.)
		return InfNaN::nan();
	x /= scale;

	if (std::isnan(x)) /* eg. original x = scale = +Inf */
		return x;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (alph == 0.) /* limit case; useful e.g. in pnchisq() */
		return (x <= 0) ? R_DT_0 : R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
	return cdfRaw(x, alph, lower_tail, log_p);
}

double Gamma::chisqQuantileAppr(double p, double nu, double g,
	bool lower_tail, bool log_p, double tol)
{
	if (std::isnan(p) || std::isnan(nu))
		return p + nu;

	constexpr auto C7 = 4.67;
	constexpr auto C8 = 6.66;
	constexpr auto C9 = 6.73;
	constexpr auto C10 = 13.32;

	double a = 0.0, ch = 0.0, p1 = 0.0;
	double p2 = 0.0, q = 0.0, t = 0.0, x = 0.0;

	if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
		return InfNaN::nan();

	if (nu <= 0) return InfNaN::nan();

	double alpha = 0.5 * nu;/* = [pq]gamma() shape */
	double c = alpha - 1;

	double t1 = log_p ? p : log(p);
	double t2 = p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p));
	double t3 = log_p ? t2 : log1p(-p);
	double t4 = lower_tail ? t1 : t3;
	double t5 = lower_tail ? t3 : t1;

	if (nu < (-1.24) * (p1 = t4)) {
		double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
		ch = exp((lgam1pa + p1) / alpha + M_LN2);
	}
	else if (nu > 0.32) {	/*  using Wilson and Hilferty estimate */
		x = Normal::quantile(p, 0, 1, lower_tail, log_p);
		p1 = 2. / (9 * nu);
		ch = nu * pow(x * sqrt(p1) + 1 - p1, 3);

		/* approximation for p tending to 1: */
		if (ch > 2.2 * nu + 6)
			ch = -2 * (t5 - c * log(0.5 * ch) + g);
	}
	else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */

		ch = 0.4;
		a = t5 + g + c * M_LN2;

		do {
			q = ch;
			p1 = 1. / (1 + ch * (C7 + ch));
			p2 = ch * (C9 + ch * (C8 + ch));
			t = -0.5 + (C7 + 2 * ch) * p1 - (C9 + ch * (C10 + 3 * ch)) / p2;
			ch -= (1 - exp(a + 0.5 * ch) * p2 * p1) / t;
		} while (fabs(q - ch) > tol * fabs(ch));
	}
	return ch;
}

double Gamma::quantile(double p, double alpha, double scale,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(alpha) || std::isnan(scale))
		return p + alpha + scale;

	constexpr auto EPS1 = 1e-2;
	constexpr auto EPS2 = 5e-7/* final precision of AS 91 */;
	constexpr auto EPS_N = 1e-15/* precision of Newton step / iterations */;
	// constexpr auto LN_EPS = -36.043653389117156 /* = log(.Machine$double.eps) iff IEEE_754 */; // unuse

	constexpr auto MAXIT = 1000/* was 20 */;

	constexpr auto pMIN = 1e-100   /* was 0.000002 = 2e-6 */;
	constexpr auto pMAX = (1-1e-14)/* was (1-1e-12) and 0.999998 = 1 - 2e-6 */;

	constexpr double i420 = 1.0 / 420.0;
	constexpr double i2520 = 1.0 / 2520.0;
	constexpr double i5040 = 1.0 / 5040.0;

	/*double a, b, c, g, ch, ch0, p1;
	double p2, q, s1, s2, s3, s4, s5, s6, t, x;*/
	int max_it_Newton = 1;

    if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if(p == 0)
			return lower_tail ? InfNaN::posinf() : 0.0;
		if(p == InfNaN::neginf())
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

	if (alpha < 0 || scale <= 0) return InfNaN::nan();

	if (alpha == 0) return 0.;

	if (alpha < 1e-10) {
		std::cout << "Warning: value of shape " << alpha
			<< " is extremely small: results may be unreliable." << std::endl;
		max_it_Newton = 7;/* may still be increased below */
	}

	double tmp = lower_tail ? (p) : (0.5 - p + 0.5);
	double p_ = log_p ? (lower_tail ? exp(p) : -expm1(p)) : tmp;
	double g = SpecialFunctions::Gamma::lgammafn(alpha);
	double ch = chisqQuantileAppr(p, 2 * alpha,  g, lower_tail, log_p, EPS1);

	double c = 0.0, s6 = 0.0, ch0 = 0.0;

	if (!std::isfinite(ch)) {
		max_it_Newton = 0; goto END;
	}
	if (ch < EPS2) {
		max_it_Newton = 20;
		goto END;
	}

	if (p_ > pMAX || p_ < pMIN) {
		max_it_Newton = 20;
		goto END;
	}

	c = alpha - 1;
	s6 = (120 + c * (346 + 127 * c)) * i5040;

	ch0 = ch;/* save initial approx. */
	for (int i = 1; i <= MAXIT; i++) {
		double q = ch;
		double p1 = 0.5 * ch;
		double p2 = p_ - cdfRaw(p1, alpha, true, false);

		if (!std::isfinite(p2) || ch <= 0)
		{
			ch = ch0; max_it_Newton = 27; goto END;
		}

		double t = p2 * exp(alpha * M_LN2 + g + p1 - c * log(ch));
		double b = t / ch;
		double a = 0.5 * t - b * c;
		double s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) * i420;
		double s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) * i2520;
		double s3 = (210 + a * (462 + a * (707 + 932 * a))) * i2520;
		double s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) * i5040;
		double s5 = (84 + 2264 * a + c * (1175 + 606 * a)) * i2520;

		ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
		if (fabs(q - ch) < EPS2 * ch)
			goto END;
		if (fabs(q - ch) > 0.1 * ch) {
			if (ch < q) ch = 0.9 * q; else ch = 1.1 * q;
		}
	}

END:
	double x = 0.5 * scale * ch;
	if (max_it_Newton) {
		/* always use log scale */
		if (!log_p) {
			p = log(p);
			log_p = true;
		}
		if (x == 0) {
			const double _1_p = 1. + 1e-7;
			const double _1_m = 1. - 1e-7;
			x = DBL_MIN;
			p_ = cdf(x, alpha, scale, lower_tail, log_p);
			if ((lower_tail && p_ > p * _1_p) ||
				(!lower_tail && p_ < p * _1_m))
				return 0.0;
		}
		else
			p_ = cdf(x, alpha, scale, lower_tail, log_p);
		if (p_ == InfNaN::neginf()) return 0;
		for (int i = 1; i <= max_it_Newton; i++) {
			double p1 = p_ - p;

			if (fabs(p1) < fabs(EPS_N * p))
				break;
			
			double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
			if ((g = pdf(x, alpha, scale, log_p)) == R_D__0) {
				break;
			}

			double t = log_p ? p1 * exp(p_ - g) : p1 / g;/* = "delta x" */
			t = lower_tail ? x - t : x + t;
			p_ = cdf(t, alpha, scale, lower_tail, log_p);
			if (fabs(p_ - p) > fabs(p1) ||
				(i > 1 && fabs(p_ - p) == fabs(p1))) {
				break;
			}
			x = t;
		}
	}
	return x;
}

double Gamma::rand(const double shape, const double scale)
{
	if (std::isnan(shape) || std::isnan(scale)) {
		return InfNaN::nan();
	}
    if (shape <= 0.0 || scale <= 0.0) {
		if(scale == 0.0 || shape == 0.0)
			return 0.0;
		return InfNaN::nan();
    }
    if(!std::isfinite(shape) || !std::isfinite(scale)) {
    	return InfNaN::posinf();
    }

    std::random_device d;   // non-deterministic random number
    std::mt19937_64 e(d()); // random engine
    std::gamma_distribution<double> u(shape, scale);
    return u(e);
}