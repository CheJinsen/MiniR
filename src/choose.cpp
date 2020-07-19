/*
 * This file is part of MiniR.
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2004-2014 The R Foundation
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

#include "specfun.h"
using namespace SpecialFunctions;

double Choose::choose(double n, double k)
{
    double r, k0 = k;
    k = nearbyint(k);

    if (std::isnan(n) || std::isnan(k)) return n + k;

    // #ifndef MATHLIB_STANDALONE
        // R_CheckStack();
    // #endif
    if (fabs(k - k0) > 1e-7) {
        std::cout << "Waring: k " << k0 << " must be integer, rounded to "
            << k << std::endl;
    }
    if (k < 30) {
        if (n - k < k && n >= 0 && isInt(n))
            k = nearbyint(n - k);
        if (k < 0) return 0.0;
        if (k == 0) return 1.0;
     
        r = n;
        for (int j = 2; j <= k; j++)
            r *= (n - j + 1) / j;
        return isInt(n) ? nearbyint(r) : r;
    }
    if (n < 0) {
        r = choose(-n + k - 1, k);
        if (isOdd(k)) r = -r;
        return r;
    }
    else if (isInt(n)) {
        n = nearbyint(n);
        if (n < k) return 0.0;
        if (n - k < 30) return choose(n, n - k);
        return nearbyint(exp(lfastchoose(n, k)));
    }
    if (n < k - 1) {
        int s_choose;
        r = lfastchoose2(n, k, &s_choose);
        return s_choose * exp(r);
    }
    return exp(lfastchoose(n, k));
}

double Choose::lchoose(double n, double k)
{
    double k0 = k;
    k = nearbyint(k);

    if (std::isnan(n) || std::isnan(k)) return n + k;

    // #ifndef MATHLIB_STANDALONE
        // R_CheckStack();
    // #endif
    if (fabs(k - k0) > 1e-7) {
        std::cout << "Warning: k " << k0 << " must be integer, rounded to "
            << k << std::endl;
    }
    if (k < 2) {
        if (k < 0) return InfNaN::neginf();
        if (k == 0) return 0.0;
        return log(fabs(n));
    }
    if (n < 0) {
        return lchoose(-n + k - 1, k);
    }
    else if (isInt(n)) {
        n = nearbyint(n);
        if (n < k) return InfNaN::neginf();
        if (n - k < 2) return lchoose(n, n - k);

        return lfastchoose(n, k);
    }
    if (n < k - 1) {
        int s = 0;
        return lfastchoose2(n, k, &s);
    }
    return lfastchoose(n, k);
}

double Choose::lfastchoose(double n, double k)
{
	return -log(n + 1.0) - Beta::lbeta(n - k + 1.0, k + 1.0);
}

double Choose::lfastchoose2(double n, double k, int* s_choose)
{
	double r = Gamma::lgammafnSign(n - k + 1.0, s_choose);
	return Gamma::lgammafn(n + 1.0) - Gamma::lgammafn(k + 1.0) - r;
}

bool Choose::isOdd(const double k)
{
	return k != 2.0 * floor(k / 2.0);
}

bool Choose::isInt(const double x)
{
	return fabs(x - nearbyint(x)) < 1e-7 * std::max(1.0, fabs(x));
}