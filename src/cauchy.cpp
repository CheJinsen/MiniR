#include "randist.h"
using namespace Randist;

double Cauchy::pdf(double x, double location, double scale, bool give_log)
{
    /* NaNs propagated correctly */
    if (std::isnan(x) || std::isnan(location) || std::isnan(scale))
        return x + location + scale;

    if (scale <= 0) return InfNaN::nan();

    double y = (x - location) / scale;
    return give_log ?
        -log(M_PI * scale * (1. + y * y)) :
        1. / (M_PI * scale * (1. + y * y));
}

double Cauchy::cdf(double x, double location, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(location) || std::isnan(scale))
        return x + location + scale;

    if (scale <= 0) return InfNaN::nan();

    x = (x - location) / scale;
    if (std::isnan(x)) return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (!std::isfinite(x)) {
        if (x < 0) return R_DT_0;
        else return R_DT_1;
    }

    if (!lower_tail)
        x = -x;

    if (fabs(x) > 1) {
        double y = atan(1 / x) / M_PI;
        double temp1 = log_p ? log1p(-y) : (0.5 - y +0.5);
        return (x > 0) ? temp1 : (log_p ? log(-y) : -y);// R_D_val(-y);
    }
    else {
        double temp2 = 0.5 + atan(x) / M_PI;
        return log_p ? log(temp2) : temp2;
    }
}

double Cauchy::quantile(double p, double location, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(location) || std::isnan(scale))
        return p + location + scale;

    if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
        return InfNaN::nan();
    if (scale <= 0 || !std::isfinite(scale)) {
        if (scale == 0) return location;
        return InfNaN::nan();
    }

    double my_INF = location + (lower_tail ? scale : -scale) * InfNaN::posinf();
    if (log_p) {
        if (p > -1) {
            if (p == 0.0)
                return my_INF;
            lower_tail = !lower_tail;
            p = -expm1(p);
        }
        else {
            p = exp(p);
        }
    }
    else {
        if (p > 0.5) {
            if (p == 1.0)
                return my_INF;
            p = 1 - p;
            lower_tail = !lower_tail;
        }
    }

    if (p == 0.5) return location; // avoid 1/Inf below
    if (p == 0.) return location + (lower_tail ? scale : -scale) * InfNaN::neginf(); // p = 1. is handled above
    return location + (lower_tail ? -scale : scale) / tanpi(p);
}

double Cauchy::tanpi(double x)
{
    if (std::isnan(x)) return x;
    if (!std::isfinite(x)) return InfNaN::nan();

    x = fmod(x, 1.);
    if (x <= -0.5) x++; else if (x > 0.5) x--;
    return (x == 0.) ? 0. : ((x == 0.5) ? InfNaN::nan() : tan(M_PI * x));
}