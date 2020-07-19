/*
 * This file is part of MiniR.
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2006 The R Core Team
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

double Uniform::pdf(double x, double a, double b, bool give_log)
{
    if (std::isnan(x) || std::isnan(a) || std::isnan(b))
        return x + a + b;

    if (b <= a) InfNaN::nan();

    if (a <= x && x <= b)
        return give_log ? -log(b - a) : 1. / (b - a);
    return give_log ? InfNaN::neginf() : 0.0;
}

double Uniform::cdf(double x, double a, double b, bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(a) || std::isnan(b))
        return x + a + b;

    if (b < a) return InfNaN::nan();
    if (!std::isfinite(a) || !std::isfinite(b)) return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (x >= b)
        return R_DT_1;
    if (x <= a)
        return R_DT_0;
    if (lower_tail) {
        double temp = (x - a) / (b - a);
        return log_p ? log(temp) : temp;
    }
    else {
        double temp = (b - x) / (b - a);
        return log_p ? log(temp) : (temp);
    }
}

double Uniform::quantile(double p, double a, double b,
    bool lower_tail, bool log_p)
{
    if (std::isnan(p) || std::isnan(a) || std::isnan(b))
        return p + a + b;

    if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
        return InfNaN::nan();

    if (!std::isfinite(a) || !std::isfinite(b)) return InfNaN::nan();
    if (b < a) return InfNaN::nan();
    if (b == a) return a;

    double temp = log_p ? (lower_tail ? exp(p) : -expm1(p))
        : (lower_tail ? p : (0.5 - p + 0.5));

    return a + temp * (b - a);
}

double Uniform::rand(const double a, const double b)
{
    if (!std::isfinite(a) || !std::isfinite(b) || b < a) {
        return InfNaN::nan();
    }

    std::random_device d;   // non-deterministic random number
    std::mt19937_64 e(d()); // random engine
    std::uniform_real_distribution<double> u(a, b);
    return u(e);
}