/*
 * This file is part of MiniR.
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2015 The R Core Team
 * Copyright (C) 2005 The R Foundation
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

double Weibull::pdf(double x, double shape, double scale, bool give_log)
{
    if (std::isnan(x) || std::isnan(shape) || std::isnan(scale))
        return x + shape + scale;
    if (shape <= 0 || scale <= 0)
        return InfNaN::nan();
    if (x < 0 || !std::isfinite(x))
        return give_log ? InfNaN::neginf() : 0.0;
    if (x == 0 && shape < 1)
        return InfNaN::posinf();

    double tmp1 = pow(x / scale, shape - 1);
    double tmp2 = tmp1 * (x / scale);
    /* These are incorrect if tmp1 == 0 */
    return  give_log ?
        -tmp2 + log(shape * tmp1 / scale) :
        shape * tmp1 * exp(-tmp2) / scale;
}

double Weibull::cdf(double x, double shape, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(shape) || std::isnan(scale))
        return x + shape + scale;
    if (shape <= 0 || scale <= 0)
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    if (x <= 0)
        return lower_tail ? R_D__0 : R_D__1;;

    x = -pow(x / scale, shape);
    const double tmp1 = x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x));
    const double tmp2 = log_p ? x : exp(x);
    return lower_tail ? (log_p ? tmp1 : -expm1(x)) : tmp2;
}

double Weibull::quantile(double p, double shape, double scale,
    bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(shape) || std::isnan(scale))
        return p + shape + scale;
    if (shape <= 0 || scale <= 0)
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

    const double tmp1 = log_p ? (p) : log(p);
    const double tmp2 = p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p));
    const double tmp3 = log_p ? tmp2 : log1p(-p);
    const double tmp4 = lower_tail ? tmp3 : tmp1;

    return scale * pow(-tmp4, 1.0 / shape);
}

double Weibull::rand(const double shape, const double scale)
{
    if (!std::isfinite(shape) || !std::isfinite(scale) ||
        shape <= 0.0 || scale <= 0.0) {
        if(scale == 0.0)
            return 0.0;
        return InfNaN::nan();
    }

    std::random_device d;   // non-deterministic random number
    std::mt19937_64 e(d()); // random engine
    std::weibull_distribution<double> u(shape, scale);
    return u(e);
}