/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 1999-2018 The R Core Team
 * Copyright (C) 2004-2018 The R Foundation
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

#include "specfun.h"

using namespace SpecialFunctions;

double Gamma::gammafn(const double x)
{
	const double gamcs[42] = {
		+0.8571195590989331421920062399942e-2,
		+0.4415381324841006757191315771652e-2,
		+0.5685043681599363378632664588789e-1,
		-0.4219835396418560501012500186624e-2,
		+0.1326808181212460220584006796352e-2,
		-0.1893024529798880432523947023886e-3,
		+0.3606925327441245256578082217225e-4,
		-0.6056761904460864218485548290365e-5,
		+0.1055829546302283344731823509093e-5,
		-0.1811967365542384048291855891166e-6,
		+0.3117724964715322277790254593169e-7,
		-0.5354219639019687140874081024347e-8,
		+0.9193275519859588946887786825940e-9,
		-0.1577941280288339761767423273953e-9,
		+0.2707980622934954543266540433089e-10,
		-0.4646818653825730144081661058933e-11,
		+0.7973350192007419656460767175359e-12,
		-0.1368078209830916025799499172309e-12,
		+0.2347319486563800657233471771688e-13,
		-0.4027432614949066932766570534699e-14,
		+0.6910051747372100912138336975257e-15,
		-0.1185584500221992907052387126192e-15,
		+0.2034148542496373955201026051932e-16,
		-0.3490054341717405849274012949108e-17,
		+0.5987993856485305567135051066026e-18,
		-0.1027378057872228074490069778431e-18,
		+0.1762702816060529824942759660748e-19,
		-0.3024320653735306260958772112042e-20,
		+0.5188914660218397839717833550506e-21,
		-0.8902770842456576692449251601066e-22,
		+0.1527474068493342602274596891306e-22,
		-0.2620731256187362900257328332799e-23,
		+0.4496464047830538670331046570666e-24,
		-0.7714712731336877911703901525333e-25,
		+0.1323635453126044036486572714666e-25,
		-0.2270999412942928816702313813333e-26,
		+0.3896418998003991449320816639999e-27,
		-0.6685198115125953327792127999999e-28,
		+0.1146998663140024384347613866666e-28,
		-0.1967938586345134677295103999999e-29,
		+0.3376448816585338090334890666666e-30,
		-0.5793070335782135784625493333333e-31
	};

	/*
	 * For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 * (xmin, xmax) are non-trivial, see ./gammalims.c
	 * xsml = exp(.01)*DBL_MIN
	 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
	*/
	const int ngam = 22;
	const double xmin = -170.5674972726612;
	const double xmax = 171.61447887182298;
	const double xsml = 2.2474362225598545e-308;
	const double dxrel = 1.490116119384765696e-8;


	if (std::isnan(x)) return x;

	if (x == 0 || (x < 0 && x == round(x))) {
		std::cerr << "Argument out of domain in gamma." << std::endl;
		return InfNaN::nan();
	}

	int n = 0;
	double y = fabs(x);
	double value = 0.0;

	if (y <= 10) {
		n = (int)x;
		if (x < 0) --n;
		y = x - n;
		--n;
		value = chebyshevEval(y * 2 - 1, gamcs, ngam) + 0.9375;
		if (n == 0)
			return value;

		if (n < 0) {
			if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
				std::cerr << "Full precision may not have been achieved in gamma." << std::endl;
			}
			if (y < xsml) {
				std::cerr << "Value out of range in gamma." << std::endl;
				if (x > 0) return InfNaN::posinf();
				else return InfNaN::neginf();
			}

			n = -n;
			for (int i = 0; i < n; i++) {
				value /= (x + i);
			}
			return value;
		}
		else {
			for (int i = 1; i <= n; i++) {
				value *= (y + i);
			}
			return value;
		}
	}
	else {
		if (x > xmax) {
			return InfNaN::posinf();
		}

		if (x < xmin) {
			return 0.0;
		}

		if (y <= 50 && y == (int)y) { /* compute (n - 1)! */
			value = 1.;
			for (int i = 2; i < y; i++) value *= i;
		}
		else { /* normal case */
			value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
				((2 * y == (int)(2 * y)) ? stirlerr(y) : lgammacor(y)));
		}
		if (x > 0)
			return value;

		if (fabs((x - (int)(x - 0.5)) / x) < dxrel) {
			std::cerr << "Full precision may not have been achieved in gamma." << std::endl;
		}

		double sinpiy = sinpi(y);
		if (sinpiy == 0) {		/* Negative integer arg - overflow */
			std::cerr << "Value out of range in gamma." << std::endl;
			return InfNaN::posinf();
		}
		return -M_PI / (y * sinpiy * value);
	}
}

double Gamma::lgammafnSign(double x, int* sgn)
{
	/*
	 * For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 * xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
	 * dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
	 */
	 
	constexpr auto xmax = 2.5327372760800758e+305;
	constexpr auto dxrel = 1.490116119384765625e-8;

	if (sgn != NULL) *sgn = 1;

	if (std::isnan(x)) return x;

	if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
		*sgn = -1;

	if (x <= 0 && x == trunc(x)) {
		return InfNaN::posinf();
	}

	double y = fabs(x);

	if (y < 1e-306) return -log(y); // denormalized range, R change
	if (y <= 10) return log(fabs(gammafn(x)));
	/*
	  ELSE  y = |x| > 10 ---------------------- */

	if (y > xmax) {
		// No warning: +Inf is the best answer
		return InfNaN::posinf();
	}

	if (x > 0) {
		if (x > 1e17)
			return(x * (log(x) - 1.0));
		else if (x > 4934720.0)
			return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
		else
			return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
	}
	/* else: x < -10; y = -x */
	double sinpiy = fabs(sinpi(y));

	if (sinpiy == 0) {
		std::cerr << " ** should NEVER happen! *** "
			<< "[gamma.cpp: Neg.int, y= " << y << std::endl;
		return InfNaN::nan();
	}

	double ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) -
		x - log(sinpiy) - lgammacor(y);

	if (fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {
		std::cerr << "Full precision may not have been "
			<< "achieved in lgamma." << std::endl;
	}

	return ans;
}

double Gamma::lgammafn(double x)
{
	return lgammafnSign(x, NULL);
}

double Gamma::lgammacor(const double x)
{
	const double algmcs[15] = {
		+0.1666389480451863247205729650822e+0,
		-0.1384948176067563840732986059135e-4,
		+0.9810825646924729426157171547487e-8,
		-0.1809129475572494194263306266719e-10,
		+0.6221098041892605227126015543416e-13,
		-0.3399615005417721944303330599666e-15,
		+0.2683181998482698748957538846666e-17,
		-0.2868042435334643284144622399999e-19,
		+0.3962837061046434803679306666666e-21,
		-0.6831888753985766870111999999999e-23,
		+0.1429227355942498147573333333333e-24,
		-0.3547598158101070547199999999999e-26,
		+0.1025680058010470912000000000000e-27,
		-0.3401102254316748799999999999999e-29,
		+0.1276642195630062933333333333333e-30
	};

	/*
	 * For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 * xbig = 2 ^ 26.5
	 * xmax = DBL_MAX / 48 =  2^1020 / 3
	 */
	const int nalgm = 5;
	const double xbig = 94906265.62425156;
	const double xmax = 3.745194030963158e306;

	if (x < 10)
		return InfNaN::nan();
	else if (x >= xmax) {
		std::cerr << "Underflow occurred in lgammacor." << std::endl;
	}
	else if (x < xbig) {
		double tmp = 10 / x;
		return chebyshevEval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
	}
	return 1 / (x * 12);
}

// sin(pi * x)  -- exact when x = k/2  for all integer k
double Gamma::sinpi(double x) {

	if (std::isnan(x)) return x;

	if (!std::isfinite(x)) return InfNaN::nan();

	x = fmod(x, 2.0);

	if (x <= -1)
		x += 2.0;
	else if (x > 1.0)
		x -= 2.0;
	if (x == 0.0 || x == 1.0) return 0.0;
	if (x == 0.5)	return  1.0;
	if (x == -0.5)	return -1.0;

	return sin(M_PI * x);
}

int Gamma::chebyshevInit(const double* dos, const int nos, const double eta)
{
    if (nos < 1)
        return 0;

    double err = 0.0;
    int i = 0;			/* just to avoid compiler warnings */
    for (int ii = 1; ii <= nos; ii++) {
        i = nos - ii;
        err += fabs(dos[i]);
        if (err > eta) {
            return i;
        }
    }
    return i;
}

double Gamma::chebyshevEval(const double x, const double* a, const int n)
{
    if (n < 1 || n > 1000) EXIT_FAILURE;
    if (x < -1.1 || x > 1.1) EXIT_FAILURE;

    double twox = x * 2;
    double b2 = 0.0;
    double b1 = 0.0;
    double b0 = 0.0;
    for (int i = 1; i <= n; i++) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

double Gamma::logspaceAdd(const double logx, const double logy)
{
	return std::max(logx, logy) + log1p(exp(-fabs(logx - logy)));
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
		return(Gamma::lgammafn(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI);
	}

	nn = n * n;
	if (n > 500) return (S0 - S1 / nn) / n;
	if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
	if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
	/* 15 < n <= 35 : */
	return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}