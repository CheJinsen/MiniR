#include "randist.h"
#include <climits>
using namespace Randist;

double Hyper::BinomialPdfRaw(const double x, const double n, const double p,
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

double Hyper::bd0(const double x, const double np)
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

double Hyper::stirlerr(const double n)
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

bool Hyper::isNegOrNonInt(double x)
{
	double isNonInt = fabs((x)- nearbyint(x)) >
		1e-7 * std::max(1.0, fabs(x));
	return x < 0.0 || isNonInt;
}

double Hyper::pdf(double x, double r, double b, double n, bool give_log)
{
	if (std::isnan(x) || std::isnan(r) || std::isnan(b) || std::isnan(n))
		return x + r + b + n;
	if (isNegOrNonInt(r) || isNegOrNonInt(b) || isNegOrNonInt(n) || n > r + b)
		return InfNaN::nan();

	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
	double R_D__1 = give_log ? 0.0 : 1.0;
	if (x < 0)
		return(R_D__0);
	
	double isNonInt = fabs((x)-nearbyint(x)) >
		1e-7 * std::max(1.0, fabs(x));
	if (isNonInt) {
		std::cout << "Warning: non-integer x = " << x << std::endl;
		return R_D__0;
	}

	x = nearbyint(x);
	r = nearbyint(r);
	b = nearbyint(b);
	n = nearbyint(n);

	if (n < x || r < x || n - x > b) return R_D__0;
	if (n == 0) return((x == 0) ? R_D__1 : R_D__0);

	double p = ((double)n) / ((double)(r + b));
	double q = ((double)(r + b - n)) / ((double)(r + b));
	double p1 = BinomialPdfRaw(x, r, p, q, give_log);
	double p2 = BinomialPdfRaw(n - x, b, p, q, give_log);
	double p3 = BinomialPdfRaw(n, r + b, p, q, give_log);

	return (give_log) ? p1 + p2 - p3 : p1 * p2 / p3;
}

double Hyper::pdhyper(double x, double NR, double NB, double n, bool log_p)
{
	/*
	 * Calculate
	 *
	 *	    phyper (x, NR, NB, n, TRUE, FALSE)
	 *   [log]  ----------------------------------
	 *	       dhyper (x, NR, NB, n, FALSE)
	 *
	 * without actually calling phyper.  This assumes that
	 *
	 *     x * (NR + NB) <= n * NR
	 *
	 */
	long double sum = 0.0;
	long double term = 1.0;

	while (x > 0 && term >= DBL_EPSILON * sum) {
		term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
		sum += term;
		x--;
	}

	double ss = (double)sum;
	return log_p ? log1p(ss) : 1 + ss;
}

double Hyper::cdf(double x, double NR, double NB, double n,
	bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(NR) || std::isnan(NB) || std::isnan(n))
		return x + NR + NB + n;

	x = floor(x + 1e-7);
	NR = nearbyint(NR);
	NB = nearbyint(NB);
	n = nearbyint(n);

	if (NR < 0 || NB < 0 || !std::isfinite(NR + NB) || n < 0 || n > NR + NB)
		return InfNaN::nan();

	if (x * (NR + NB) > n * NR) {
		// Swap tails.
		double oldNB = NB;
		NB = NR;
		NR = oldNB;
		x = n - x - 1;
		lower_tail = !lower_tail;
	}

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (x < 0)
		return R_DT_0;
	if (x >= NR || x >= n)
		return R_DT_1;

	double d = Hyper::pdf(x, NR, NB, n, log_p);
	double pd = pdhyper(x, NR, NB, n, log_p);

	double x2 = d + pd;
	double tmp2 = d * pd;
	double t1 = x2 > -M_LN2 ? log(-expm1(x2)) : log1p(-exp(x2));
	double t2 = lower_tail ? x2 : t1;
	double t3 = lower_tail ? tmp2 : (0.5 - tmp2 + 0.5);

	return log_p ? t2 : t3;
}

double Hyper::quantile(double p, double NR, double NB, double n,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(NR) || std::isnan(NB) || std::isnan(n))
		return p + NR + NB + n;

	if (!std::isfinite(p) || !std::isfinite(NR) || !std::isfinite(NB) || !std::isfinite(n))
		return InfNaN::nan();

	NR = nearbyint(NR);
	NB = nearbyint(NB);
	double N = NR + NB;
	n = nearbyint(n);
	if (NR < 0 || NB < 0 || n < 0 || n > N)
		return InfNaN::nan();

	double xstart = std::max(0.0, n - NB);
	double xend = std::min(n, NR);

	//R_Q_P01_boundaries(p, xstart, xend);
	if (log_p) {
		if (p > 0)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? xend : xstart;
		if (p == InfNaN::neginf())
			return lower_tail ? xstart : xend;
	}
	else {
		if (p < 0 || p > 1)
			return InfNaN::nan();
		if (p == 0)
			return lower_tail ? xstart : xend;
		if (p == 1)
			return lower_tail ? xend : xstart;
	}

	double xr = xstart;
	double xb = n - xr;
	double sum = 0.0;
	bool small_N = (N < 1000);
	double term = SpecialFunctions::Choose::lfastchoose(NR, xr) +
		SpecialFunctions::Choose::lfastchoose(NB, xb) -
		SpecialFunctions::Choose::lfastchoose(N, n);
	if (small_N)
		term = exp(term);
	NR -= xr;
	NB -= xb;

	double t1 = lower_tail ? p : (0.5 - p + 0.5);
	double t2 = log_p ? (lower_tail ? exp(p) : -expm1(p)) : t1;
	if (!lower_tail || log_p) {
		p = t2;
	}
	p *= 1 - 1000 * DBL_EPSILON;
	sum = small_N ? term : exp(term);

	while (sum < p && xr < xend) {
		xr++;
		NB++;
		if (small_N)
			term *= (NR / xr) * (xb / NB);
		else
			term += log((NR / xr) * (xb / NB));
		sum += small_N ? term : exp(term);
		xb--;
		NR--;
	}
	return xr;
}

double Hyper::afc(int i)
{
    // If (i > 7), use Stirling's approximation, otherwise use table lookup.
    const double al[8] = {
		0.0,/*ln(0!)=ln(1)*/
		0.0,/*ln(1!)=ln(1)*/
		0.69314718055994530941723212145817,/*ln(2) */
		1.79175946922805500081247735838070,/*ln(6) */
		3.17805383034794561964694160129705,/*ln(24)*/
		4.78749174278204599424770093452324,
		6.57925121201010099506017829290394,
		8.52516136106541430016553103634712
    };

    if (i < 0) {
		std::cout << "Warning: Hyper::rand(): afc(i), i = " << i
			<< " < 0 -- SHOULD NOT HAPPEN!" << std::endl;
		return -1;
    }
    if (i <= 7){
		return al[i];
    }
    
    double di = i, i2 = di*di;
    return (di + 0.5) * log(di) - di + M_LN_SQRT_2PI +
		(0.0833333333333333 - 0.00277777777777778 / i2) / di;
}

//     rhyper(NR, NB, n) -- NR 'red', NB 'blue', n drawn, how many are 'red'
double Hyper::rand(double nn1in, double nn2in, double kkin)
{
    int ks = -1, n1s = -1, n2s = -1;
    int k = 0, n1 = 0, n2 = 0;

    // II :
    double w = 0.0;
    // III:
    double a = 0.0, d = 0.0, s = 0.0, xl = 0.0, xr = 0.0;
    double kl = 0.0, kr = 0.0, lamdl = 0.0, lamdr = 0.0;
    double p1 = 0.0, p2 = 0.0, p3 = 0.0;

    if(!std::isfinite(nn1in) || !std::isfinite(nn2in) || !std::isfinite(kkin)) {
		return InfNaN::nan();
    }

    nn1in = nearbyint(nn1in);
    nn2in = nearbyint(nn2in);
    kkin  = nearbyint(kkin);

    if (nn1in < 0 || nn2in < 0 || kkin < 0 || kkin > nn1in + nn2in) {
    	return InfNaN::nan();
    }
    if (nn1in >= INT_MAX || nn2in >= INT_MAX || kkin >= INT_MAX) {
		if(kkin == 1.) {
		    return Binomial::rand(kkin, nn1in / (nn1in + nn2in));
		}
		return quantile(Uniform::rand(), nn1in, nn2in, kkin, false, false);
    }

    int nn1 = (int)nn1in;
    int nn2 = (int)nn2in;
    int kk  = (int)kkin;
    bool setup1 = false, setup2 = false;

    if (nn1 != n1s || nn2 != n2s) { // n1 | n2 is changed: setup all
		setup1 = true;	setup2 = true;
    }
    else if (kk != ks) { // n1 & n2 are unchanged: setup 'k' only
		setup1 = false;	setup2 = true;
    }
    else { // all three unchanged ==> no setup
		setup1 = false;	setup2 = false;
    }

    double N = 0.0;
    if (setup1) { // n1 & n2
		n1s = nn1; n2s = nn2; // save
		N = nn1 + (double)nn2; // avoid int overflow
		if (nn1 <= nn2) {
		    n1 = nn1; n2 = nn2;
		}
		else { // nn2 < nn1
		    n1 = nn2; n2 = nn1;
		}
    }

    if (setup2) { // k
		ks = kk; // save
		if ((double)kk + kk >= N) { // this could overflow
		    k = (int)(N - kk);
		}
		else {
		    k = kk;
		}
    }

	int m = 0, minjx = 0, maxjx = 0;
    if (setup1 || setup2) {
		m = (int) ((k + 1.) * (n1 + 1.) / (N + 2.)); // m := floor(adjusted mean E[.])
		minjx = std::max(0.0, 1.0 * k - n2);
		maxjx = std::min(1.0 * n1, 1.0 * k);
    }

    int ix = 0;
    if (minjx == maxjx) {
		ix = maxjx;
		goto L_finis; // return appropriate variate

    }
    else if (m - minjx < 10) { // II: (Scaled) algorithm HIN (inverse transformation) ----
		constexpr double scale = 1e25; // scaling factor against (early) underflow
		constexpr double con = 57.5646273248511421;

		if (setup1 || setup2) {
		    double lw = 0.0; // log(w);  w = exp(lw) * scale = exp(lw + log(scale)) = exp(lw + con)
		    if (k < n2) {
				lw = afc(n2) + afc(n1 + n2 - k) - afc(n2 - k) - afc(n1 + n2);
		    }
		    else {
				lw = afc(n1) + afc(     k     ) - afc(k - n2) - afc(n1 + n2);
		    }
		    w = exp(lw + con);
		}
		double p = 0.0, u = 0.0;

	    L10:
		p = w;
		ix = minjx;
		u = Uniform::rand() * scale;

		while (u > p) {
		    u -= p;
		    p *= ((double) n1 - ix) * (k - ix);
		    ix++;
		    p = p / ix / (n2 - k + ix);

		    if (ix > maxjx)
				goto L10;
		}
    }
    else {

		double u = 0.0, v = 0.0;

		if (setup1 || setup2) {
		    s = sqrt((N - k) * k * n1 * n2 / (N - 1) / N / N);

		    /* remark: d is defined in reference without int. */
		    /* the truncation centers the cell boundaries at 0.5 */

		    d = (int) (1.5 * s) + .5;
		    xl = m - d + .5;
		    xr = m + d + .5;
		    a = afc(m) + afc(n1 - m) + afc(k - m) + afc(n2 - k + m);
		    kl = exp(a - afc((int)(xl)) - afc((int)(n1 - xl))
			     - afc((int)(k - xl))
			     - afc((int)(n2 - k + xl)));
		    kr = exp(a - afc((int) (xr - 1))
			     - afc((int)(n1 - xr + 1))
			     - afc((int)(k - xr + 1))
			     - afc((int)(n2 - k + xr - 1)));
		    lamdl = -log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
		    lamdr = -log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
		    p1 = d + d;
		    p2 = p1 + kl / lamdl;
		    p3 = p2 + kr / lamdr;
		}

		int n_uv = 0;
	    
	    L30:
		u = Uniform::rand() * p3;
		v = Uniform::rand();
		n_uv++;

		if(n_uv >= 10000) {
		    // REprintf("rhyper(*, n1=%d, n2=%d, k=%d): branch III: giving up after %d rejections\n",
			   //   nn1, nn2, kk, n_uv);
		    return InfNaN::nan();
	    }

		if (u < p1) {		/* rectangular region */
		    ix = (int) (xl + u);
		}
		else if (u <= p2) {	/* left tail */
		    ix = (int) (xl + log(v) / lamdl);
		    if (ix < minjx)
			goto L30;
		    v = v * (u - p1) * lamdl;
		}
		else {		/* right tail */
		    ix = (int) (xr - log(v) / lamdr);
		    if (ix > maxjx)
			goto L30;
		    v = v * (u - p2) * lamdr;
		}

		/* acceptance/rejection test */
		bool reject = true;

		if (m < 100 || ix <= 50) {
		    double f = 1.0;
		    if (m < ix) {
				for (int i = m + 1; i <= ix; i++)
				    f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
		    }
		    else if (m > ix) {
				for (int i = ix + 1; i <= m; i++)
				    f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
			}
		    if (v <= f) {
				reject = false;
		    }
		}
		else {

		    constexpr double deltal = 0.0078;
		    constexpr double deltau = 0.0034;

		    double de = 0.0, dr = 0.0, ds = 0.0, dt = 0.0;

		    /* squeeze using upper and lower bounds */
		    double y = ix;
		    double y1 = y + 1.0;
		    double ym = y - m;
		    double yn = n1 - y + 1.0;
		    double yk = k - y + 1.0;
		    double nk = n2 - k + y1;
		    double r = -ym / y1;
		    s = ym / yn;
		    double t = ym / yk;
		    double e = -ym / nk;
		    double g = yn * yk / (y1 * nk) - 1.0;
		    double dg = 1.0;

		    if (g < 0.0)
				dg = 1.0 + g;

		    double gu = g * (1.0 + g * (-0.5 + g / 3.0));
		    double gl = gu - .25 * (g * g * g * g) / dg;
		    double xm = m + 0.5;
		    double xn = n1 - m + 0.5;
		    double xk = k - m + 0.5;
		    double nm = n2 - k + xm;
		    double ub = y * gu - m * gl + deltau
				+ xm * r * (1. + r * (-0.5 + r / 3.0))
				+ xn * s * (1. + s * (-0.5 + s / 3.0))
				+ xk * t * (1. + t * (-0.5 + t / 3.0))
				+ nm * e * (1. + e * (-0.5 + e / 3.0));
		    /* test against upper bound */
		    double alv = log(v);
		    if (alv > ub) {
				reject = true;
		    }
		    else {
				dr = xm * (r * r * r * r);
				if (r < 0.0)
				    dr /= (1.0 + r);

				ds = xn * (s * s * s * s);
				if (s < 0.0)
				    ds /= (1.0 + s);

				dt = xk * (t * t * t * t);
				if (t < 0.0)
				    dt /= (1.0 + t);

				de = nm * (e * e * e * e);
				if (e < 0.0)
				    de /= (1.0 + e);
				if (alv < ub - 0.25 * (dr + ds + dt + de)
				    + (y + m) * (gl - gu) - deltal) {
				    reject = false;
				}
				else {
				    if (alv <= (a - afc(ix) - afc(n1 - ix)
						- afc(k - ix) - afc(n2 - k + ix))) {
						reject = false;
				    }
				    else {
						reject = true;
				    }
				}
		    }
		} // else
		if (reject)
		    goto L30;

    } // end{branch III}


L_finis:  // return appropriate variate

    if ((double)kk + kk >= N) {
		if (nn1 > nn2) {
		    ix = kk - nn2 + ix;
		} else {
		    ix = nn1 - ix;
		}
    }
    else if (nn1 > nn2) {
		ix = kk - ix;
    }

    return ix;
}