#include "randist.h"
using namespace Randist;

double Fdist::BinomialPdfRaw(const double x, const double n, const double p,
	const double q, bool log_p)
{
	double lf = 0.0, lc = 0.0;
	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;

	if (p == 0) return((x == 0) ? R_D__1 : R_D__0);
	if (q == 0) return((x == n) ? R_D__1 : R_D__0);

	if (x == 0) {
		if (n == 0) return R_D__1;
		lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q);
		return log_p ? lc : exp(lc);
	}
	if (x == n) {
		lc = (q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p);
		return log_p ? lc : exp(lc);
	}
	if (x < 0 || x > n) return(R_D__0);

	lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);

	lf = M_LN_2PI + log(x) + log1p(-x / n);

	return log_p ? (lc - 0.5 * lf) : exp(lc - 0.5 * lf);
}

double Fdist::bd0(const double x, const double np)
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

double Fdist::stirlerr(const double n)
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
		dens = BinomialPdfRaw((m - 2) / 2, (m + n - 2) / 2, p, q, give_log);
	}
	else {
		f = m * m * q / (2 * p * (m + n));
		dens = BinomialPdfRaw(m / 2, (m + n) / 2, p, q, give_log);
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