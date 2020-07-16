#include "randist.h"

using namespace Randist;

double Normal::pdf(double x, double mu, double sigma, bool give_log)
{
    if (std::isnan(x) || std::isnan(mu) || std::isnan(sigma))
        return x + mu + sigma;

    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;

    if (sigma < 0) return InfNaN::nan();
    if (!std::isfinite(sigma)) return R_D__0;
    if (!std::isfinite(x) && mu == x) return InfNaN::nan();/* x-mu is NaN */
    if (sigma == 0)
        return (x == mu) ? InfNaN::posinf() : R_D__0;
    x = (x - mu) / sigma;

    if (!std::isfinite(x)) return R_D__0;

    x = fabs(x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
        return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));

    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

	double temp = 1.0 * DBL_MIN_EXP + 1.0 - 1.0 * DBL_MANT_DIG;
    if (x > sqrt(-2.0 * M_LN2 * temp))
        return 0.0;

    double x1 = ldexp(nearbyint(ldexp(x, 16)), -16);
    double x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
        (exp(-0.5 * x1 * x1) * exp((-0.5 * x2 - x1) * x2));
}

double Normal::cdf(double x, double mu, double sigma,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(mu) || std::isnan(sigma))
		return x + mu + sigma;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (!std::isfinite(x) && mu == x) return InfNaN::nan();/* x-mu is NaN */
	if (sigma <= 0) {
		if (sigma < 0) return InfNaN::nan();
		return (x < mu) ? R_DT_0 : R_DT_1;
	}
	double p = (x - mu) / sigma;
	if (!std::isfinite(p))
		return (x < mu) ? R_DT_0 : R_DT_1;
	x = p;

	double cp = 0.0;
	cdfBoth(x, p, cp, (lower_tail ? 0 : 1), log_p);

	return(lower_tail ? p : cp);
}

void Normal::cdfBoth(double x, double& cum, double& ccum,
	bool i_tail, bool log_p)
{
	const double a[5] = {
		2.2352520354606839287,
		161.02823106855587881,
		1067.6894854603709582,
		18154.981253343561249,
		0.065682337918207449113
	};
	const double b[4] = {
		47.20258190468824187,
		976.09855173777669322,
		10260.932208618978205,
		45507.789335026729956
	};
	const double c[9] = {
		0.39894151208813466764,
		8.8831497943883759412,
		93.506656132177855979,
		597.27027639480026226,
		2494.5375852903726711,
		6848.1904505362823326,
		11602.651437647350124,
		9842.7148383839780218,
		1.0765576773720192317e-8
	};
	const double d[8] = {
		22.266688044328115691,
		235.38790178262499861,
		1519.377599407554805,
		6485.558298266760755,
		18615.571640885098091,
		34900.952721145977266,
		38912.003286093271411,
		19685.429676859990727
	};
	const double p[6] = {
		0.21589853405795699,
		0.1274011611602473639,
		0.022235277870649807,
		0.001421619193227893466,
		2.9112874951168792e-5,
		0.02307344176494017303
	};
	const double q[5] = {
		1.28426009614491121,
		0.468238212480865118,
		0.0659881378689285515,
		0.00378239633202758244,
		7.29751555083966205e-5
	};

	double xden = 0.0, xnum = 0.0, temp = 0.0;
	double del = 0.0, xsq = 0.0;

	if (std::isnan(x)) { cum = ccum = x; return; }

	double eps = DBL_EPSILON * 0.5;

	/* i_tail in {0,1,2} =^= {lower, upper, both} */
	bool lower = (i_tail != 1);
	bool upper = (i_tail != 0);

	double y = fabs(x);
	if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
		if (y > eps) {
			xsq = x * x;
			xnum = a[4] * xsq;
			xden = xsq;
			for (int i = 0; i < 3; ++i) {
				xnum = (xnum + a[i]) * xsq;
				xden = (xden + b[i]) * xsq;
			}
		}
		else xnum = xden = 0.0;

		temp = x * (xnum + a[3]) / (xden + b[3]);
		if (lower)  cum = 0.5 + temp;
		if (upper) ccum = 0.5 - temp;
		if (log_p) {
			if (lower)  cum = log(cum);
			if (upper) ccum = log(ccum);
		}
	}
	else if (y <= M_SQRT_32) {

		/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

		xnum = c[8] * y;
		xden = y;
		for (int i = 0; i < 7; ++i) {
			xnum = (xnum + c[i]) * y;
			xden = (xden + d[i]) * y;
		}
		temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * 16) / 16;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    ccum = 1.0 - cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = cum; if(lower) cum = ccum; ccum = temp;	\
	}

		do_del(y);
		swap_tail;
	}
	else if ((log_p && y < 1e170)
		|| (lower && -37.5193 < x && x < 8.2924)
		|| (upper && -8.2924 < x && x < 37.5193)) {

		xsq = 1.0 / (x * x);
		xnum = p[5] * xsq;
		xden = xsq;
		for (int i = 0; i < 4; ++i) {
			xnum = (xnum + p[i]) * xsq;
			xden = (xden + q[i]) * xsq;
		}
		temp = xsq * (xnum + p[4]) / (xden + q[4]);
		temp = (M_1_SQRT_2PI - temp) / y;

		do_del(x);
		swap_tail;
	}
	else {
		if (x > 0) {
			cum = log_p ? 0.0 : 1.0;
			ccum = log_p ? InfNaN::neginf() : 0.0;
		}
		else {
			cum = log_p ? InfNaN::neginf() : 0.0;
			ccum = log_p ? 0.0 : 1.0;
		}
	}
	return;
}

double Normal::quantile(double p, double mu, double sigma,
	bool lower_tail, bool log_p)
{
	double r = 0.0, val = 0.0;

	if (std::isnan(p) || std::isnan(mu) || std::isnan(sigma))
		return p + mu + sigma;

	//R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
	if (log_p) {
		if (p > 0)	return InfNaN::nan();
		if (p == 0) return lower_tail ? InfNaN::posinf() : InfNaN::neginf();
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

	if (sigma < 0)	return InfNaN::nan();
	if (sigma == 0)	return mu;

	double p_ = log_p ? (lower_tail ? exp(p) : -expm1(p))
		: (lower_tail ? (p) : (0.5 - p + 0.5));
	double q = p_ - 0.5;

	/*
	 * -- use AS 241 ---
	 * double ppnd16_(double *p, long *ifault)
	 * ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
	 */
	if (fabs(q) <= 0.425) {
		r = 0.180625 - q * q;
		val =
			q * (((((((r * 2509.0809287301226727 +
				33430.575583588128105) * r + 67265.770927008700853) * r +
				45921.953931549871457) * r + 13731.693765509461125) * r +
				1971.5909503065514427) * r + 133.14166789178437745) * r +
				3.387132872796366608)
			/ (((((((r * 5226.495278852854561 +
				28729.085735721942674) * r + 39307.89580009271061) * r +
				21213.794301586595867) * r + 5394.1960214247511077) * r +
				687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
	}
	else {
		if (q > 0)
			r = log_p ? (lower_tail ? -expm1(p) : exp(p))
			: (lower_tail ? (0.5 - p + 0.5) : p);
		else
			r = p_;

		r = sqrt(-((log_p &&
			((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : log(r)));

		if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
			r += -1.6;
			val = (((((((r * 7.7454501427834140764e-4 +
				0.0227238449892691845833) * r + 0.24178072517745061177) *
				r + 1.27045825245236838258) * r +
				3.64784832476320460504) * r + 5.7694972214606914055) *
				r + 4.6303378461565452959) * r +
				1.42343711074968357734)
				/ (((((((r *
					1.05075007164441684324e-9 + 5.475938084995344946e-4) *
					r + 0.0151986665636164571966) * r +
					0.14810397642748007459) * r + 0.68976733498510000455) *
					r + 1.6763848301838038494) * r +
					2.05319162663775882187) * r + 1.);
		}
		else { /* very close to  0 or 1 */
			r += -5.0;
			val = (((((((r * 2.01033439929228813265e-7 +
				2.71155556874348757815e-5) * r +
				0.0012426609473880784386) * r + 0.026532189526576123093) *
				r + 0.29656057182850489123) * r +
				1.7848265399172913358) * r + 5.4637849111641143699) *
				r + 6.6579046435011037772)
				/ (((((((r *
					2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
					r + 1.8463183175100546818e-5) * r +
					7.868691311456132591e-4) * r + 0.0148753612908506148525)
					* r + 0.13692988092273580531) * r +
					0.59983220655588793769) * r + 1.);
		}

		if (q < 0.0)
			val = -val;
	}
	return mu + sigma * val;
}

unsigned long int Normal::seed = time(0);

void Normal::setSeed(const unsigned long int s)
{
	seed = s;
}

double Normal::rand(const double mu, const double sigma)
{
	if (std::isnan(mu) || !std::isfinite(sigma) || sigma < 0.0)  {
		return InfNaN::nan();
	}

	std::random_device d;	// non-deterministic random number
	std::mt19937_64 e(d());	// random engine
	std::normal_distribution<double> u(mu, sigma);

    if (sigma == 0.0 || !std::isfinite(mu)) {
		return mu;
    }
    else {
    	return u(e);
    }
}