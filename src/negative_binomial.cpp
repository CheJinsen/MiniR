#include "randist.h"
using namespace Randist;

double NegBinomial::BinomialPdfRaw(const double x, const double n, const double p,
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

double NegBinomial::bd0(const double x, const double np)
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

double NegBinomial::stirlerr(const double n)
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

double NegBinomial::poissonPdfRaw(double x, double lambda, bool give_log)
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

double NegBinomial::pdf(double x, double size, double prob, bool give_log)
{
    if (std::isnan(x) || std::isnan(size) || std::isnan(prob))
        return x + size + prob;
    if (prob <= 0 || prob > 1 || size < 0)
        return InfNaN::nan();

    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
    double R_D__1 = give_log ? 0.0 : 1.0;

    bool isNonInt = fabs((x)-nearbyint(x)) > 1e-7 * std::max(1.0, fabs(x));
    if (isNonInt) {
        std::cout << "Warning: non-integer x = " << x << std::endl;
        return give_log ? InfNaN::neginf() : 0.0;
    }

    if (x < 0 || !std::isfinite(x)) return R_D__0;
    /* limiting case as size approaches zero is point mass at zero */
    if (x == 0 && size == 0) return R_D__1;
    x = nearbyint(x);
    if (!std::isfinite(size)) size = DBL_MAX;

    double ans = BinomialPdfRaw(size, x + size, prob, 1 - prob, give_log);
    double p = ((double)size) / (size + x);
    return((give_log) ? log(p) + ans : p * ans);
}

double NegBinomial::pdf_mu(double x, double size, double mu, bool give_log)
{
    if (std::isnan(x) || std::isnan(size) || std::isnan(mu)) {
        return x + size + mu;
    }
	if (mu < 0 || size < 0) {
		return InfNaN::nan();
	}

	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
	double R_D__1 = give_log ? 0.0 : 1.0;

	double isNonInt = fabs((x)-nearbyint(x)) >
		1e-7 * std::max(1.0, fabs(x));
	if (isNonInt) {
		std::cout << "Warning: non-integer x = " << x << std::endl;
		return R_D__0;
	}
    if (x < 0 || !std::isfinite(x)) return R_D__0;
    if (x == 0 && size == 0) return R_D__1;

    x = nearbyint(x);
    if (!std::isfinite(size)) { // limit case: Poisson
        return(poissonPdfRaw(x, mu, give_log));
    }

	if (x == 0) {
		double temp = size * (size < mu ? log(size / (size + mu))
			: log1p(-mu / (size + mu)));
		return give_log ? temp : exp(temp);
	}

    if (x < 1e-10 * size) {
        double p = (size < mu ? log(size / (1 + size / mu))
			: log(mu / (1 + mu / size)));
		double temp = x * p - mu - lgamma(x + 1) +
			log1p(x * (x - 1) / (2 * size));
        return give_log ? temp : exp(temp);
    }
    else {
		double p = ((double)size) / (size + x);
        double ans = BinomialPdfRaw(size, x + size, size / (size + mu),
			mu / (size + mu), give_log);
        return (give_log) ? log(p) + ans : p * ans;
    }
}

double NegBinomial::cdf(double x, double size, double prob,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(size) || std::isnan(prob))
		return x + size + prob;
	if (!std::isfinite(size) || !std::isfinite(prob))	return InfNaN::nan();
	if (size < 0 || prob <= 0 || prob > 1)	return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (size == 0)
		return (x >= 0) ? R_DT_1 : R_DT_0;

	if (x < 0) return R_DT_0;
	if (!std::isfinite(x)) return R_DT_1;
	x = floor(x + 1e-7);
	return Beta::cdf(prob, size, x + 1, lower_tail, log_p);
}

double NegBinomial::cdf_mu(double x, double size, double mu,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(size) || std::isnan(mu))
		return x + size + mu;
	if (!std::isfinite(mu))	return InfNaN::nan();
	if (size < 0 || mu < 0)	return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (size == 0)
		return (x >= 0) ? R_DT_1 : R_DT_0;

	if (x < 0) return R_DT_0;
	if (!std::isfinite(x)) return R_DT_1;
	if (!std::isfinite(size))
		return Poisson::cdf(x, mu, lower_tail, log_p);

	x = floor(x + 1e-7);
	/* return
	 * pbeta(pr, size, x + 1, lower_tail, log_p);  pr = size/(size + mu), 1-pr = mu/(size+mu)
	 *
	 *= pbeta_raw(pr, size, x + 1, lower_tail, log_p)
	 *            x.  pin   qin
	 *=  bratio (pin,  qin, x., 1-x., &w, &wc, &ierr, log_p),  and return w or wc ..
	 *=  bratio (size, x+1, pr, 1-pr, &w, &wc, &ierr, log_p) */
	{
		int ierr = 0;
		double w = 0.0, wc = 0.0;
		Toms::bratio(size, x + 1, size / (size + mu), mu / (size + mu), &w, &wc, &ierr, log_p);
		if (ierr)
			std::cout << "Warning: cdf_mu() -> bratio() gave error code " << ierr << std::endl;
		return lower_tail ? w : wc;
	}
}

double NegBinomial::doSearch(double y, double* z, double p, double n,
	double pr, double incr)
{
	if (*z >= p) {	/* search to the left */
		for (;;) {
			if (y == 0 ||
				(*z = cdf(y - incr, n, pr, true, false)) < p)
				return y;
			y = std::max(0.0, y - incr);
		}
	}
	else {		/* search to the right */
		for (;;) {
			y = y + incr;
			if ((*z = cdf(y, n, pr, true, false)) >= p)
				return y;
		}
	}
}

double NegBinomial::quantile(double p, double size, double prob,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(size) || std::isnan(prob))
		return p + size + prob;

	if (prob == 0 && size == 0) return 0;

	if (prob <= 0 || prob > 1 || size < 0) return InfNaN::nan();

	if (prob == 1 || size == 0) return 0;

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

	double Q = 1.0 / prob;
	double P = (1.0 - prob) * Q;
	double mu = size * P;
	double sigma = sqrt(size * P * Q);
	double gamma = (Q + P) / sigma;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	const double temp1 = lower_tail ? (p) : (0.5 - (p)+0.5);
	if (!lower_tail || log_p) {
		p = log_p ? (lower_tail ? exp(p) : -expm1(p)) : temp1;
		if (p == R_DT_0) return 0;
		if (p == R_DT_1) return InfNaN::posinf();
	}

	if (p + 1.01 * DBL_EPSILON >= 1.0) return InfNaN::posinf();

	double z = Normal::quantile(p, 0.0, 1.0, true, false);
	double y = nearbyint(mu + sigma * (z + gamma * (z * z - 1) / 6));

	z = cdf(y, size, prob, true, false);
	p *= 1 - 64 * DBL_EPSILON;

	/* If the C-F value is not too large a simple search is OK */
	if (y < 1e5)
		return doSearch(y, &z, p, size, prob, 1);

	{
		double incr = floor(y * 0.001), oldincr;
		do {
			oldincr = incr;
			y = doSearch(y, &z, p, size, prob, incr);
			incr = std::max(1.0, floor(incr / 100.0));
		} while (oldincr > 1 && incr > y * 1e-15);
		return y;
	}
}

double NegBinomial::quantile_mu(double p, double size, double mu,
	bool lower_tail, bool log_p)
{
	if (size == InfNaN::posinf())
		return(Poisson::quantile(p, mu, lower_tail, log_p));
	return quantile(p, size, size / (size + mu), lower_tail, log_p);
}