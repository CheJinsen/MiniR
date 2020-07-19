/*
 * This file is part of MiniR.
 * Copyright (C) 1999-2014  The R Core Team
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

std::vector<double> Signrank::w;
int Signrank::allocated_n = 0;

void Signrank::wInitMaybe(int n)
{
    int u = n * (n + 1) / 2;
    int c = (u / 2);

    w.clear();

    if(w.empty()) {
        w.resize(c + 1);
	    allocated_n = n;
    }
}

double Signrank::csignrank(int k, int n)
{
// #ifndef MATHLIB_STANDALONE
//     R_CheckUserInterrupt();
// #endif

    int u = n * (n + 1) / 2;
    int c = (u / 2);

    if (k < 0 || k > u)
	   return 0;
    if (k > c)
	   k = u - k;

    if (n == 1)
        return 1.0;
    if (w[0] == 1.0)
        return w[k];

    w[0] = w[1] = 1.0;
    for(int j = 2; j < n+1; ++j) {
        int end = std::min(j * (j + 1.0) / 2.0, 1.0 * c);

    	for(int i = end; i >= j; --i)
    	    w[i] += w[i-j];
    }

    return w[k];
}

double Signrank::dtvalue(const double x, const bool lower_tail, const bool log_p)
{
    double t1 = log_p ? log1p(-x) : (0.5 - x + 0.5);
    double t2 = log_p ? log(x) : x;
    return lower_tail ? t2 : t1;
}

double Signrank::pdf(double x, double n, bool give_log)
{
    if (std::isnan(x) || std::isnan(n))
        return x + n;

    n = nearbyint(n);
    if (n <= 0)
	   return InfNaN::nan();

    if (fabs(x - nearbyint(x)) > 1e-7)
	   return give_log ? InfNaN::neginf() : 0.0;
    x = nearbyint(x);
    if ((x < 0) || (x > (n * (n + 1) / 2)))
	   return give_log ? InfNaN::neginf() : 0.0;

    int nn = (int)n;
    wInitMaybe(nn);

    double temp = log(csignrank((int)x, nn)) - n * M_LN2;
    double d = give_log ? temp : exp(temp);

    return d;
}

double Signrank::cdf(double x, double n, bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(n))
        return x + n;

    if (!std::isfinite(n))
        return InfNaN::nan();

    n = nearbyint(n);
    if (n <= 0)
        return InfNaN::nan();

    x = nearbyint(x + 1e-7);
    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (x < 0.0)
	   return R_DT_0;
    if (x >= n * (n + 1) / 2)
	   return R_DT_1;

    int nn = (int) n;
    wInitMaybe(nn);
    double f = exp(- n * M_LN2);
    double p = 0.0;

    if (x <= (n * (n + 1) / 4)) {
    	for (int i = 0; i <= x; i++)
    	    p += csignrank(i, nn) * f;
    }
    else {
    	x = n * (n + 1) / 2 - x;
    	for (int i = 0; i < x; i++)
    	    p += csignrank(i, nn) * f;
    	lower_tail = !lower_tail;
    }

    return dtvalue(p, lower_tail, log_p);
}

double Signrank::quantile(double x, double n, bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(n))
	   return x + n;

    if (!std::isfinite(x) || !std::isfinite(n))
	   return InfNaN::nan();

    // R_Q_P01_check(x);
    if ((log_p  && x > 0) || (!log_p && (x < 0 || x > 1)))
        return InfNaN::nan();

    n = nearbyint(n);
    if (n <= 0)
	   return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (x == R_DT_0)
	   return 0.0;
    if (x == R_DT_1)
	   return n * (n + 1.0) / 2.0;

    double temp = lower_tail ? x : (0.5 - x + 0.5);
    if(log_p || !lower_tail)
	   x = log_p ? (lower_tail ? exp(x) : - expm1(x)) : temp;

    int nn = (int) n;
    wInitMaybe(nn);
    double f = exp(- n * M_LN2);
    double p = 0.0;
    int q = 0;

    if (x <= 0.5) {
    	x = x - 10 * DBL_EPSILON;
    	for (;;) {
    	    p += csignrank(q, nn) * f;
    	    if (p >= x)
    	       	break;
    	    q++;
    	}
    }
    else {
    	x = 1 - x + 10 * DBL_EPSILON;
    	for (;;) {
    	    p += csignrank(q, nn) * f;
    	    if (p > x) {
        		q = (int)(n * (n + 1) / 2 - q);
        		break;
    	    }
    	    q++;
    	}
    }

    return q;
}

double Signrank::rand(double n)
{
    if (std::isnan(n))
        return n;

    n = nearbyint(n);
    if (n < 0)
        return InfNaN::nan();

    if (n == 0)
        return 0;
    
    double r = 0.0;
    int k = (int)n;
    for (int i = 0; i < k; ) {
        r += (++i) * floor(Uniform::rand() + 0.5);
    }
    return r;
}