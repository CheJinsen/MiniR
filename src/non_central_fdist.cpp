#include "randist.h"
using namespace Randist;

double NonCentralFdist::pdf(double x, double df1, double df2,
    double ncp, bool give_log)
{
    if (std::isnan(x) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
        return x + df2 + df1 + ncp;
    if (df1 <= 0. || df2 <= 0. || ncp < 0)
        return InfNaN::nan();
    if (x < 0.0)
        return give_log ? InfNaN::neginf() : 0.0;
    if (!std::isfinite(ncp))
        return InfNaN::nan();

    if (!std::isfinite(df1) && !std::isfinite(df2)) {
        if (x == 1.0)
            return InfNaN::posinf();
        else
            return give_log ? InfNaN::neginf() : 0.0;
    }
    if (!std::isfinite(df2))
        return df1 * NonCentralChisq::pdf(x * df1, df1, ncp, give_log);
    
    double f = 0.0;
    double z = 0.0;
    if (df1 > 1e14 && ncp < 1e7) {
        f = 1 + ncp / df1;
        z = Gamma::cdf(1.0 / x / f, df2 / 2, 2.0 / df2, give_log);
        return give_log ? z - 2 * log(x) - log(f) : z / (x * x) / f;
    }

    double y = (df1 / df2) * x;
    z = NonCentralBeta::pdf(y / (1 + y), df1 / 2.0, df2 / 2.0, ncp, give_log);
    return  give_log ?
        z + log(df1) - log(df2) - 2 * log1p(y) :
        z * (df1 / df2) / (1 + y) / (1 + y);
}

double NonCentralFdist::cdf(double x, double df1, double df2, double ncp,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
        return x + df2 + df1 + ncp;
    if (df1 <= 0. || df2 <= 0. || ncp < 0)
        return InfNaN::nan();
    if (!std::isfinite(ncp))
        return InfNaN::nan();
    if (!std::isfinite(df1) && !std::isfinite(df2)) /* both +Inf */
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
    if (x <= 0.0) return R_DT_0;
    if (x >= InfNaN::posinf()) return R_DT_1;

    if (df2 > 1e8)
        return NonCentralChisq::cdf(x * df1, df1, ncp, lower_tail, log_p);

    double y = (df1 / df2) * x;
    return NonCentralBetaCdf2(y / (1. + y), 1. / (1. + y), df1 / 2., df2 / 2.,
        ncp, lower_tail, log_p);
}

long double NonCentralFdist::NonCentralBetaCdfRaw(double x, double o_x, double a,
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

double NonCentralFdist::NonCentralBetaCdf2(double x, double o_x, double a, double b, double ncp,
	bool lower_tail, bool log_p)
{
	long double ans = NonCentralBetaCdfRaw(x, o_x, a, b, ncp);

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

double NonCentralFdist::quantile(double p, double df1, double df2, double ncp,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
		return p + df1 + df2 + ncp;
	if (df1 <= 0. || df2 <= 0. || ncp < 0)
		return InfNaN::nan();
	if (!std::isfinite(ncp))
		return InfNaN::nan();
	if (!std::isfinite(df1) && !std::isfinite(df2))
		return InfNaN::nan();

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

	if (df2 > 1e8) /* avoid problems with +Inf and loss of accuracy */
		return NonCentralChisq::quantile(p, df1, ncp, lower_tail, log_p) / df1;

	double y = NonCentralBeta::quantile(p, df1 / 2.0, df2 / 2.0, ncp, lower_tail, log_p);
	return y / (1 - y) * (df2 / df1);
}