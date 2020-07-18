#include "randist.h"
using namespace Randist;

double Exp::pdf(double x, double scale, bool give_log)
{
    /* NaNs propagated correctly */
    if (std::isnan(x) || std::isnan(scale)) return x + scale;

    if (scale <= 0.0) return InfNaN::nan();

    if (x < 0.0)
        return give_log ? InfNaN::neginf() : 0.0;
    return (give_log ?
        (-x / scale) - log(scale) :
        exp(-x / scale) / scale);
}

double Exp::cdf(double x, double scale, bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(scale))
        return x + scale;
    if (scale < 0) return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    if (x <= 0.0)
        return R_DT_0;

    /* same as weibull( shape = 1): */
    x = -(x / scale);
    double temp1 = x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x));
    return lower_tail
        ? (log_p ? temp1 : -expm1(x))
        : (log_p ? (x) : exp(x));
}

double Exp::quantile(double p, double scale, bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(scale))
        return p + scale;

    if (scale < 0) return InfNaN::nan();

    if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    if (p == R_DT_0)
        return 0;

    double temp1 = log_p ? (p) : log(p);
    double temp2 = p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p));
    double temp3 = log_p ? temp2 : log1p(-p);
    return -scale * (lower_tail ? temp3 : temp1);
}

double Exp::rand(const double lambda)
{
    if (!std::isfinite(lambda) || lambda <= 0.0) {
        if(lambda == 0.0)
            return 0.0;
        return InfNaN::nan();
    }

    std::random_device d;   // non-deterministic random number
    std::mt19937_64 e(d()); // random engine
    std::exponential_distribution<double> u(lambda);
    return u(e);
}