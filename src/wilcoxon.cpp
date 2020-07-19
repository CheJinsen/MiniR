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

int Wilcox::allocated_m = 0;
int Wilcox::allocated_n = 0;
std::vector<std::vector<std::vector<double>>> Wilcox::w;


void Wilcox::wInitMaybe(int m, int n)
{
    if (m > n) {
        int temp = n; n = m; m = temp;
    }

    w.clear();

    if (w.empty()) {
        m = std::max(m, 50);
        n = std::max(n, 50);
        w.resize(m + 1);

        for (int i = 0; i <= m; i++) {
            w[i].resize(n + 1);
        }
        allocated_m = m; allocated_n = n;
    }
}


// This counts the number of choices with statistic = k
double Wilcox::cwilcox(int k, int m, int n)
{
//#ifndef MATHLIB_STANDALONE
//    R_CheckUserInterrupt();
//#endif

    int u = m * n;
    if (k < 0 || k > u)
        return 0;
    int c = (int)(u / 2);
    if (k > c)
        k = u - k;

    int i = 0, j = 0;
    if (m < n) {
        i = m; j = n;
    }
    else {
        i = n; j = m;
    }

    if (j == 0)
        return k == 0;

    if (j > 0 && k < j) return cwilcox(k, i, k);


    if (w[i][j].size() == 0) {
        w[i][j].resize(c + 1);

        for (int l = 0; l <= c; l++)
            w[i][j][l] = -1;
    }
    if (w[i][j][k] < 0) {
        if (j == 0)
            w[i][j][k] = (k == 0);
        else
            w[i][j][k] = cwilcox(k - j, i - 1, j) + cwilcox(k, i, j - 1);

    }
    return w[i][j][k];
}

double Wilcox::pdf(double x, double m, double n, bool give_log)
{
    if (std::isnan(x) || std::isnan(m) || std::isnan(n))
        return x + m + n;

    m = nearbyint(m);
    n = nearbyint(n);
    if (m <= 0 || n <= 0)
        return InfNaN::nan();

    if (fabs(x - nearbyint(x)) > 1e-7)
        return give_log ? InfNaN::neginf() : 0.0;
    x = nearbyint(x);
    if ((x < 0) || (x > m * n))
        return give_log ? InfNaN::neginf() : 0.0;

    int mm = (int)m, nn = (int)n, xx = (int)x;
    wInitMaybe(mm, nn);
    double d = give_log ?
        log(cwilcox(xx, mm, nn)) - SpecialFunctions::Choose::lchoose(m + n, n) :
        cwilcox(xx, mm, nn) / SpecialFunctions::Choose::choose(m + n, n);

    return d;
}

double Wilcox::dtvalue(const double x, const bool lower_tail, const bool log_p)
{
    double t1 = log_p ? log1p(-x) : (0.5 - x + 0.5);
    double t2 = log_p ? log(x) : x;
    return lower_tail ? t2 : t1;
}

double Wilcox::cdf(double q, double m, double n,
    bool lower_tail, bool log_p)
{
    if (std::isnan(q) || std::isnan(m) || std::isnan(n))
        return q + m + n;

    if (!std::isfinite(m) || !std::isfinite(n))
        return InfNaN::nan();
    m = nearbyint(m);
    n = nearbyint(n);
    if (m <= 0 || n <= 0)
        return InfNaN::nan();

    q = floor(q + 1e-7);
    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (q < 0.0)
        return R_DT_0;
    if (q >= m * n)
        return R_DT_1;

    int mm = (int)m, nn = (int)n;
    wInitMaybe(mm, nn);
    double c = SpecialFunctions::Choose::choose(m + n, n);
    double p = 0.0;

    if (q <= (m * n / 2)) {
        for (int i = 0; i <= q; i++)
            p += cwilcox(i, mm, nn) / c;
    }
    else {
        q = m * n - q;
        for (int i = 0; i < q; i++)
            p += cwilcox(i, mm, nn) / c;
        lower_tail = !lower_tail;
    }

    return dtvalue(p, lower_tail, log_p);
}

double Wilcox::quantile(double x, double m, double n,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(m) || std::isnan(n))
        return x + m + n;
    if (!std::isfinite(x) || !std::isfinite(m) || !std::isfinite(n))
        return InfNaN::nan();
    
    //R_Q_P01_check(x);
    if ((log_p && x > 0) || (!log_p && (x < 0 || x > 1)))
        return InfNaN::nan();

    m = nearbyint(m);
    n = nearbyint(n);
    if (m <= 0 || n <= 0)
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

    if (x == R_DT_0)
        return 0;
    if (x == R_DT_1)
        return m * n;

    double temp = lower_tail ? x : (0.5 - x + 0.5);
    if (log_p || !lower_tail)
        x = log_p ? (lower_tail ? exp(x) : -expm1(x)) : temp;

    int mm = (int)m, nn = (int)n;
    wInitMaybe(mm, nn);
    double c = SpecialFunctions::Choose::choose(m + n, n);
    double p = 0.0;
    int q = 0;
    if (x <= 0.5) {
        x = x - 10 * DBL_EPSILON;
        for (;;) {
            p += cwilcox(q, mm, nn) / c;
            if (p >= x)
                break;
            q++;
        }
    }
    else {
        x = 1 - x + 10 * DBL_EPSILON;
        for (;;) {
            p += cwilcox(q, mm, nn) / c;
            if (p > x) {
                q = (int)(m * n - q);
                break;
            }
            q++;
        }
    }
    return q;
}

//generate a random non-negative integer < 2 ^ bits in 16 bit chunks
double Wilcox::rbits(const int bits)
{
    long long v = 0;
    for (int n = 0; n <= bits; n += 16) {
        int v1 = (int)floor(Uniform::rand() * 65536);
        v = 65536 * v + v1;
    }
    const long long one64 = 1L;
    // mask out the bits in the result that are not needed
    return (double)(v & ((one64 << bits) - 1));
}

double Wilcox::uniformIndex(const double dn)
{
    // rejection sampling from integers below the next larger power of two
    if (dn <= 0) {
        return 0.0;
    }
    int bits = (int)ceil(log2(dn));
    double dv = 0.0;
    do { dv = rbits(bits); } while (dn <= dv);
    return dv;
}

double Wilcox::rand(double m, double n)
{
    if (std::isnan(m) || std::isnan(n))
        return m + n;

    m = nearbyint(m);
    n = nearbyint(n);
    if ((m < 0) || (n < 0))
        return InfNaN::nan();

    if ((m == 0) || (n == 0))
        return 0;

    double r = 0.0;
    int k = (int) (m + n);
    std::vector<int> x;

    for (int i = 0; i < k; i++) {
        x.push_back(i);
    }

    int j = 0;
    for (int i = 0; i < n; i++) {
        j = (int)uniformIndex(k);
        r += x[j];
        x[j] = x[--k];
    }
    
    return r - n * (n - 1) / 2;
}