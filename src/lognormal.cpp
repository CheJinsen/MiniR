#include "randist.h"
using namespace Randist;

double Lognormal::pdf(double x, double meanlog, double sdlog, bool give_log)
{
    double y = 0.0;

    if (std::isnan(x) || std::isnan(meanlog) || std::isnan(sdlog))
        return x + meanlog + sdlog;

    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;
    if (sdlog < 0) InfNaN::nan();
    if (!std::isfinite(x) && log(x) == meanlog)
        return InfNaN::nan();/* log(x) - meanlog is NaN */
    if (sdlog == 0)
        return (log(x) == meanlog) ? InfNaN::posinf() : R_D__0;
    if (x <= 0) return R_D__0;

    y = (log(x) - meanlog) / sdlog;
    return give_log ?
        -(M_LN_SQRT_2PI + 0.5 * y * y + log(x * sdlog)) :
        M_1_SQRT_2PI * exp(-0.5 * y * y) / (x * sdlog);
}

double Lognormal::cdf(double x, double meanlog, double sdlog,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(meanlog) || std::isnan(sdlog))
        return x + meanlog + sdlog;

    if (sdlog < 0) return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;

    if (x > 0)
        return Normal::cdf(log(x), meanlog, sdlog, lower_tail, log_p);
    return R_DT_0;
}

double Lognormal::quantile(double p, double meanlog, double sdlog,
    bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(meanlog) || std::isnan(sdlog))
        return p + meanlog + sdlog;

    if (log_p) {
        if (p > 0) {
            return InfNaN::nan();
        }
	    if(p == 0) {
	        return lower_tail ? InfNaN::posinf() : 0.0;
        }
        if (p == InfNaN::neginf()) {
            return lower_tail ? 0.0 : InfNaN::posinf();
        }
    }							
    else {
        if (p < 0 || p > 1) {
            return InfNaN::nan();
        }
	    if(p == 0) {
	        return lower_tail ? 0.0 : InfNaN::posinf();
        }
	    if(p == 1) {
	        return lower_tail ? InfNaN::posinf() : 0.0;
        }
    }

    return exp(Normal::quantile(p, meanlog, sdlog, lower_tail, log_p));
}