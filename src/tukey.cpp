#include "randist.h"
using namespace Randist;

double Tukey::wprob(double w, double rr, double cc)
{
	constexpr auto nleg = 12;
	constexpr auto ihalf = 6;
	constexpr double C1 = -30.0;
	constexpr double C2 = -50.0;
	constexpr double C3 = 60.0;
	constexpr double bb = 8.0;
	constexpr double wlar = 3.0;
	constexpr double wincr1 = 2.0;
	constexpr double wincr2 = 3.0;

	const double xleg[ihalf] = {
		0.981560634246719250690549090149,
		0.904117256370474856678465866119,
		0.769902674194304687036893833213,
		0.587317954286617447296702418941,
		0.367831498998180193752691536644,
		0.125233408511468915472441369464
	};
	const double aleg[ihalf] = {
		0.047175336386511827194615961485,
		0.106939325995318430960254718194,
		0.160078328543346226334652529543,
		0.203167426723065921749064455810,
		0.233492536538354808760849898925,
		0.249147045813402785000562436043
	};

	double qsqz = w * 0.5;
	if (qsqz >= bb)
		return 1.0;

	double pr_w = 2 * Normal::cdf(qsqz, 0.0, 1.0, true, false) - 1.0;

	if (pr_w >= exp(C2 / cc))
		pr_w = pow(pr_w, cc);
	else
		pr_w = 0.0;

	double wincr = 0.0;
	if (w > wlar)
		wincr = wincr1;
	else
		wincr = wincr2;

	long double blb = qsqz;
	double binc = (bb - qsqz) / wincr;
	long double bub = blb + binc;
	long double einsum = 0.0;

	/* integrate over each interval */
	double wi = 0.0;
	double cc1 = cc - 1.0;
	for (wi = 1; wi <= wincr; wi++) {
		long double elsum = 0.0;
		double a = (double)(0.5 * (bub + blb));

		/* legendre quadrature with order = nleg */

		double xx = 0.0;
		double b = (double)(0.5 * (bub - blb));
		for (int jj = 1; jj <= nleg; jj++) {
			int j = 0;
			if (ihalf < jj) {
				j = (nleg - jj) + 1;
				xx = xleg[j - 1];
			}
			else {
				j = jj;
				xx = -xleg[j - 1];
			}
			double c = b * xx;
			double ac = a + c;
			double qexpo = ac * ac;
			if (qexpo > C3)
				break;

			double pplus = 2 * Normal::cdf(ac, 0., 1., 1, 0);
			double pminus = 2 * Normal::cdf(ac, w, 1., 1, 0);
			double rinsum = (pplus * 0.5) - (pminus * 0.5);
			if (rinsum >= exp(C1 / cc1)) {
				rinsum = (aleg[j - 1] * exp(-(0.5 * qexpo))) * pow(rinsum, cc1);
				elsum += rinsum;
			}
		}
		elsum *= (((2.0 * b) * cc) * M_1_SQRT_2PI);
		einsum += elsum;
		blb = bub;
		bub += binc;
	}

	pr_w += (double)einsum;
	if (pr_w <= exp(C1 / rr))
		return 0.;

	pr_w = pow(pr_w, rr);
	if (pr_w >= 1.0)
		return 1.0;
	return pr_w;
}

double Tukey::dtvalue(const double x, const bool lower_tail, const bool log_p)
{
	double t1 = log_p ? log1p(-x) : (0.5 - x + 0.5);
	double t2 = log_p ? log(x) : x;
	return lower_tail ? t2 : t1;
}

double Tukey::cdf(double q, double rr, double cc, double df,
	bool lower_tail, bool log_p)
{
	constexpr auto nlegq = 16;
	constexpr auto ihalfq = 8;
	constexpr double eps1 = -30.0;
	constexpr double eps2 = 1.0e-14;
	constexpr double dhaf = 100.0;
	constexpr double dquar = 800.0;
	constexpr double deigh = 5000.0;
	constexpr double dlarg = 25000.0;
	constexpr double ulen1 = 1.0;
	constexpr double ulen2 = 0.5;
	constexpr double ulen3 = 0.25;
	constexpr double ulen4 = 0.125;

	const double xlegq[ihalfq] = {
		0.989400934991649932596154173450,
		0.944575023073232576077988415535,
		0.865631202387831743880467897712,
		0.755404408355003033895101194847,
		0.617876244402643748446671764049,
		0.458016777657227386342419442984,
		0.281603550779258913230460501460,
		0.950125098376374401853193354250e-1
	};
	const double alegq[ihalfq] = {
		0.271524594117540948517805724560e-1,
		0.622535239386478928628438369944e-1,
		0.951585116824927848099251076022e-1,
		0.124628971255533872052476282192,
		0.149595988816576732081501730547,
		0.169156519395002538189312079030,
		0.182603415044923588866763667969,
		0.189450610455068496285396723208
	};

	if (std::isnan(q) || std::isnan(rr) || std::isnan(cc) || std::isnan(df))
		return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (q <= 0)
		return R_DT_0;

	if (df < 2 || rr < 1 || cc < 2) return InfNaN::nan();

	if (!std::isfinite(q))
		return R_DT_1;

	if (df > dlarg)
		return dtvalue(wprob(q, rr, cc), lower_tail, log_p);

	double f2 = df * 0.5;
	/* lgammafn(u) = log(gamma(u)) */
	double f2lf = ((f2 * log(df)) - (df * M_LN2)) - SpecialFunctions::Gamma::lgammafn(f2);
	double f21 = f2 - 1.0;
	double ulen = 0.0;
	double ff4 = df * 0.25;

	if (df <= dhaf)			ulen = ulen1;
	else if (df <= dquar)	ulen = ulen2;
	else if (df <= deigh)	ulen = ulen3;
	else					ulen = ulen4;

	f2lf += log(ulen);

	double ans = 0.0;
	double otsum = 0.0;
	for (int i = 1; i <= 50; i++) {
		otsum = 0.0;
		int j = 0;
		double t1 = 0.0;
		double twa1 = (2.0 * i - 1.0) * ulen;

		for (int jj = 1; jj <= nlegq; jj++) {
			if (ihalfq < jj) {
				j = jj - ihalfq - 1;
				t1 = (f2lf + (f21 * log(twa1 + (xlegq[j] * ulen))))
					- (((xlegq[j] * ulen) + twa1) * ff4);
			}
			else {
				j = jj - 1;
				t1 = (f2lf + (f21 * log(twa1 - (xlegq[j] * ulen))))
					+ (((xlegq[j] * ulen) - twa1) * ff4);

			}

			double qsqz = 0.0;
			if (t1 >= eps1) {
				if (ihalfq < jj) {
					qsqz = q * sqrt(((xlegq[j] * ulen) + twa1) * 0.5);
				}
				else {
					qsqz = q * sqrt(((-(xlegq[j] * ulen)) + twa1) * 0.5);
				}

				double wprb = wprob(qsqz, rr, cc);
				double rotsum = (wprb * alegq[j]) * exp(t1);
				otsum += rotsum;
			}
		}

		if (i * ulen >= 1.0 && otsum <= eps2)
			break;
		ans += otsum;
	}

	if (otsum > eps2) { // not converged
		std::cout << "Full precision may not have been "
			<< "achieved in Tukey::cdf()" << std::endl;
	}
	if (ans > 1.)
		ans = 1.;
	return dtvalue(ans, lower_tail, log_p);
}

double Tukey::qinv(double p, double c, double v)
{
	constexpr double p0 = 0.322232421088;
	constexpr double q0 = 0.993484626060e-01;
	constexpr double p1 = -1.0;
	constexpr double q1 = 0.588581570495;
	constexpr double p2 = -0.342242088547;
	constexpr double q2 = 0.531103462366;
	constexpr double p3 = -0.204231210125;
	constexpr double q3 = 0.103537752850;
	constexpr double p4 = -0.453642210148e-04;
	constexpr double q4 = 0.38560700634e-02;
	constexpr double c1 = 0.8832;
	constexpr double c2 = 0.2368;
	constexpr double c3 = 1.214;
	constexpr double c4 = 1.208;
	constexpr double c5 = 1.4142;
	constexpr double vmax = 120.0;

	double ps = 0.5 - 0.5 * p;
	double yi = sqrt(log(1.0 / (ps * ps)));
	double t = yi + ((((yi * p4 + p3) * yi + p2) * yi + p1) * yi + p0)
		/ ((((yi * q4 + q3) * yi + q2) * yi + q1) * yi + q0);

	if (v < vmax) t += (t * t * t + t) / v / 4.0;
	double q = c1 - c2 * t;
	if (v < vmax) q += -c3 / v + c4 * t / v;
	return t * (q * log(c - 1.0) + c5);
}

double Tukey::quantile(double p, double rr, double cc, double df,
	bool lower_tail, bool log_p)
{
	const double eps = 0.0001;
	const int maxiter = 50;

	if (std::isnan(p) || std::isnan(rr) || std::isnan(cc) || std::isnan(df)) {
		std::cout << "Warngin: Argument out of domain in "
			<< "Tukey::quantile()." << std::endl;
		return p + rr + cc + df;
	}
	if (df < 2 || rr < 1 || cc < 2) return InfNaN::nan();

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

	double tmp = lower_tail ? (p) : (0.5 - (p)+0.5);
	p = log_p ? (lower_tail ? exp(p) : -expm1(p)) : tmp;
	double x0 = qinv(p, cc, df);
	double valx0 = cdf(x0, rr, cc, df, true, false) - p;
	double x1 = 0.0;
	if (valx0 > 0.0)
		x1 = std::max(0.0, x0 - 1.0);
	else
		x1 = x0 + 1.0;

	double valx1 = cdf(x1, rr, cc, df, true, false) - p;
	double ans = 0.0;
	for (int iter = 1; iter < maxiter; iter++) {
		ans = x1 - ((valx1 * (x1 - x0)) / (valx1 - valx0));
		valx0 = valx1;
		x0 = x1;
		if (ans < 0.0) {
			ans = 0.0;
			valx1 = -p;
		}

		valx1 = cdf(ans, rr, cc, df, true, false) - p;
		x1 = ans;
		double xabs = fabs(x1 - x0);
		if (xabs < eps)
			return ans;
	}

	std::cout << "Warning: convergence failed in "
		<< "Tukey::quantile()." << std::endl;
	return ans;
}