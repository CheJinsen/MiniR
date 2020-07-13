#include "randist.h"
using namespace Randist;

double Chisq::pdf(double x, double df, bool give_log)
{
    return Gamma::pdf(x, df / 2.0, 2.0, give_log);
}

double Chisq::cdf(double x, double df, bool lower_tail, bool log_p)
{
    return Gamma::cdf(x, df / 2.0, 2.0, lower_tail, log_p);
}

double Chisq::quantile(double p, double df, bool lower_tail, bool log_p)
{
    return Gamma::quantile(p, 0.5 * df, 2.0, lower_tail, log_p);
}