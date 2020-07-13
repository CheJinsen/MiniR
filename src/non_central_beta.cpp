#include "randist.h"
using namespace Randist;

double NonCentralBeta::bd0(const double x, const double np)
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

double NonCentralBeta::stirlerr(const double n)
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
	return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

double NonCentralBeta::poissonPdfRaw(double x, double lambda, bool give_log)
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

double NonCentralBeta::pdf(double x, double a, double b, double ncp, bool give_log)
{
    constexpr double eps = 1.e-15;

    if (std::isnan(x) || std::isnan(a) || std::isnan(b) || std::isnan(ncp))
        return x + a + b + ncp;
    if (ncp < 0 || a <= 0 || b <= 0)
        return InfNaN::nan();
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(ncp))
        return InfNaN::nan();

    if (x < 0 || x > 1)
        return give_log ? InfNaN::neginf() : 0.0;
    if (ncp == 0)
        return Beta::pdf(x, a, b, give_log);

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
    long double p_k = poissonPdfRaw(kMax, ncp2, true);
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

long double NonCentralBeta::cdfRaw(double x, double o_x, double a,
	double b, double ncp)
{
	constexpr double errmax = 1.0e-9;
	constexpr int itrmax = 10000;

	if (ncp < 0.0 || a <= 0.0 || b <= 0.0)
		return InfNaN::nan();
	if (x < 0.0 || o_x > 1.0 || (x == 0.0 && o_x == 1.0))
		return 0.0;
	if (x > 1.0 || o_x < 0.0 || (x == 1.0 && o_x == 0.0))
		return 1.0;

	double c = ncp / 2.0;
	double x0 = floor(std::max(c - 7.0 * sqrt(c), 0.0));
	double a0 = a + x0;
	double lbeta = SpecialFunctions::Gamma::lgammafn(a0) +
		SpecialFunctions::Gamma::lgammafn(b) -
		SpecialFunctions::Gamma::lgammafn(a0 + b);

	int ierr = 0;
	double temp = 0.0;
	double tmp_c = 0.0;
	Toms::bratio(a0, b, x, o_x, &temp, &tmp_c, &ierr, false);

	long double q = 0.0;
	long double gx = exp(a0 * log(x) + b * (x < .5 ? log1p(-x) : log(o_x))
		- lbeta - log(a0));
	if (a0 > a)
		q = exp(-c + x0 * log(c) - SpecialFunctions::Gamma::lgammafn(x0 + 1.0));
	else
		q = exp(-c);

	long double sumq = 1. - q;
	long double ans = q * temp;
	long double ax = q * temp;
	double j = floor(x0);
	double errbd = 0.0;
	do {
		j++;
		temp -= (double)gx;
		gx *= x * (a + b + j - 1.) / (a + j);
		q *= c / j;
		sumq -= q;
		ax = temp * q;
		ans += ax;
		errbd = (double)((temp - gx) * sumq);
	} while (errbd > errmax && j < itrmax + x0);

	if (errbd > errmax)
		std::cout << "Warning: Full precision may not have "
			<< "been achieved in cdf." << std::endl;
	if (j >= itrmax + x0)
		std::cout << "Warning: Convergence failed in cdf." << std::endl;

	return ans;
}

double NonCentralBeta::cdf2(double x, double o_x, double a, double b, double ncp,
	bool lower_tail, bool log_p)
{
	long double ans = cdfRaw(x, o_x, a, b, ncp);

	if (lower_tail) {
		return (double)(log_p ? logl(ans) : ans);
	}
	else {
		if (ans > 1. - 1e-10)
			std::cout << "Warning: Full precision may not have "
				<< "been achieved in cdf." << std::endl;
		if (ans > 1.0) ans = 1.0;
		return (double)(log_p ? log1pl(-ans) : (1.0 - ans));
	}
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
	return cdf2(x, 1 - x, a, b, ncp, lower_tail, log_p);
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