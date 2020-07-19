/*
 * This file is part of MiniR.
 *
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2019 The R Core Team
 * Copyright (C) 2005-2016 The R Foundation
 * Copyright (C) 2020 Jinsen Che
 *
 * MiniR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MiniR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MiniR. If not, see <https://www.gnu.org/licenses/>.
 */

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

double Lognormal::rand(const double meanlog, const double sdlog)
{
    if(std::isnan(meanlog) || !std::isfinite(sdlog) || sdlog < 0.0) {
        return InfNaN::nan();
    }

    std::random_device d;   // non-deterministic random number
    std::mt19937_64 e(d()); // random engine
    std::lognormal_distribution<double> u(meanlog, sdlog);
    return u(e);
}