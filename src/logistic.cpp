/*
 * This file is part of MiniR.
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2008 The R Core Team
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

double Logistic::pdf(double x, double location, double scale, bool give_log)
{
    if (std::isnan(x) || std::isnan(location) || std::isnan(scale))
        return x + location + scale;
    if (scale <= 0.0)
        return InfNaN::nan();

    x = fabs((x - location) / scale);
    double e = exp(-x);
    double f = 1.0 + e;
    return give_log ? -(x + log(scale * f * f)) : e / (scale * f * f);
}

double Logistic::log1pexp(const double x) {
    if (x <= 18.0)  return log1p(exp(x));
    if (x > 33.3)   return x;
    return x + exp(-x);
}

double Logistic::cdf(double x, double location, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(location) || std::isnan(scale))
        return x + location + scale;
    if (scale <= 0.0)
        return InfNaN::nan();

    x = (x - location) / scale;
    if (std::isnan(x))
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
    if (!std::isfinite(x)) {
        if (x > 0) return R_DT_1;
        return R_DT_0;
    }

    if (log_p) {
        // log(1 / (1 + exp( +- x ))) = -log(1 + exp( +- x))
        return -log1pexp(lower_tail ? -x : x);
    }
    else {
        return 1 / (1 + exp(lower_tail ? -x : x));
    }
}

double Logistic::quantile(double p, double location, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(location) || std::isnan(scale))
        return p + location + scale;

    //R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
    if (log_p) {
        if (p > 0)
            return InfNaN::nan();
        if (p == 0)
            return lower_tail ? InfNaN::posinf() : InfNaN::neginf();
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

    if (scale < 0.0) return InfNaN::nan();
    if (scale == 0.0) return location;

    const double tmp = (p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p)));
    if (log_p) {
        if (lower_tail)
            p = p - tmp;
        else
            p = tmp - p;
    }
    else {
        p = log(lower_tail ? (p / (1.0 - p)) : ((1.0 - p) / p));
    }

    return location + scale * p;
}

double Logistic::rand(const double location, const double scale)
{
    if (std::isnan(location) || !std::isfinite(scale)) {
        return InfNaN::nan();
    }

    if (scale == 0.0 || !std::isfinite(location)) {
        return location;
    }
    else {
        double u = Uniform::rand();
        return location + scale * log(u / (1. - u));
    }
}