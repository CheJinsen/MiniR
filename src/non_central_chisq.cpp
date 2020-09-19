/*
 * This file is part of MiniR.
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 1998 Ross Ihaka
 * Copyright (C) 2000-2015 The R Core Team
 * Copyright (C) 2004-2016 The R Foundation
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

double NonCentralChisq::pdf(double x, double df, double ncp, bool give_log)
{
    const double eps = 5e-15;
    double i = 0.0, q = 0.0, mid = 0.0, dfmid = 0.0;
    long double sum = 0.0, term = 0.0;
    double R_D__0 = give_log ? InfNaN::neginf() : 0.0;

    if (std::isnan(x) || std::isnan(df) || std::isnan(ncp))
        return x + df + ncp;


    if (!std::isfinite(df) || !std::isfinite(ncp) || ncp < 0 || df < 0)
        return InfNaN::nan();

    if (x < 0) return R_D__0;
    if (x == 0 && df < 2.)
        return InfNaN::posinf();
    if (ncp == 0)
        return (df > 0) ? Chisq::pdf(x, df, give_log) : R_D__0;
    if (x == InfNaN::posinf()) return R_D__0;

    double ncp2 = 0.5 * ncp;
    double imax = ceil((-(2 + df) + sqrt((2 - df) * (2 - df) + 4 * ncp * x)) / 4);
    if (imax < 0) imax = 0.0;
    if (std::isfinite(imax)) {
        dfmid = df + 2 * imax;
        mid = Base::poissonPdfRaw(imax, ncp2, false) * Chisq::pdf(x, dfmid, false);
    }
    else {
        mid = 0;
    }

    if (mid == 0) {
        if (give_log || ncp > 1000.) {
            double nl = df + ncp, ic = nl / (nl + ncp);
            return Chisq::pdf(x * ic, nl * ic, give_log);
        }
        else
            return R_D__0;
    }

    sum = mid;
    term = mid; df = dfmid; i = imax;
    double x2 = x * ncp2;
    do {
        i++;
        q = x2 / i / df;
        df += 2;
        term *= q;
        sum += term;
    } while (q >= 1 || term * q > (1 - q) * eps || term > 1e-10 * sum);

    term = mid; df = dfmid; i = imax;
    while (i != 0) {
        df -= 2;
        q = i * df / x2;
        i--;
        term *= q;
        sum += term;
        if (q < 1 && term * q <= (1 - q) * eps) break;
    }
    return give_log ? log((double)sum) : (double)sum;
}

double NonCentralChisq::cdf(double x, double df, double ncp,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(df) || std::isnan(ncp))
        return x + df + ncp;
    if (!std::isfinite(df) || !std::isfinite(ncp))
        return InfNaN::nan();

    if (df < 0.0 || ncp < 0.0) return InfNaN::nan();

    double ans = cdfRaw(x, df, ncp, 1e-12, 8 * DBL_EPSILON,
        1000000, lower_tail, log_p);
    double R_D__1 = log_p ? 0.0 : 1.0;

    if (x <= 0. || x == InfNaN::posinf())
        return ans; // because it's perfect

    if (ncp >= 80) {
        if (lower_tail) {
            ans = std::min(ans, R_D__1);
        }
        else {
            if (ans < (log_p ? (-10. * M_LN10) : 1e-10))
                std::cout << "Full precision may not have been achieved "
                    << "in Non-central cdf()" << std::endl;
            if (!log_p && ans < 0.) ans = 0.;
        }
    }

    if (!log_p || ans < -1e-8)
        return ans;
    else {
        ans = cdfRaw(x, df, ncp, 1e-12, 8 * DBL_EPSILON, 1000000, !lower_tail, false);
        return log1p(-ans);
    }
}

double NonCentralChisq::cdfRaw(double x, double f, double theta,
    double errmax, double reltol, int itrmax,
    bool lower_tail, bool log_p)
{
    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
    const double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;

    if (x <= 0.0) {
        if (x == 0. && f == 0.) {
            double L = -0.5 * theta;
            double t1 = log_p ? L : exp(L);
            double t2 = L > -M_LN2 ? log(-expm1(L)) : log1p(-exp(L));
            return lower_tail ? t1 : (log_p ? t2 : -expm1(L));
        }
        return R_DT_0;
    }
    if (!std::isfinite(x))	return R_DT_1;

//    /* This is principally for use from qnchisq */
//#ifndef MATHLIB_STANDALONE
//    R_CheckUserInterrupt();
//#endif

    if (theta < 80) {
        long double ans = 0.0;
        if (lower_tail && f > 0. &&
            log(x) < M_LN2 + 2 / f * (lgamma(f / 2. + 1) + _dbl_min_exp)) {
            double lambda = 0.5 * theta;
            double sum, sum2, pr = -lambda;
            sum = sum2 = InfNaN::neginf();
            /* we need to renormalize here: the result could be very close to 1 */
            for (int i = 0; i < 110; pr += log(lambda) - log(++i)) {
                sum2 = logspaceAdd(sum2, pr);
                sum = logspaceAdd(sum, pr + Chisq::cdf(x, f + 2.0 * i, lower_tail, true));
                if (sum2 >= -1e-15)
                   break;
            }
            ans = sum - sum2;
            return (double)(log_p ? ans : expl(ans));
        }
        else {
            long double lambda = 0.5 * theta; // < 40
            long double sum = 0, sum2 = 0, pr = expl(-lambda);

            for (int i = 0; i < 110; pr *= lambda / ++i) {
                sum2 += pr;
                sum += pr * Chisq::cdf(x, f + 2.0 * i, lower_tail, false);
                if (sum2 >= 1 - 1e-15) break;
            }
            ans = sum / sum2;
            return (double)(log_p ? logl(ans) : ans);
        }
    } // if(theta < 80)

    double u = 0.0;
    long double lu = 0.0;
    double l_lam = 0.0;
    double lam = 0.5 * theta;
    bool lamSml = (-lam < _dbl_min_exp);
    if (lamSml) {
        lu = -lam;
        l_lam = log(lam);
    }
    else {
        u = exp(-lam);
    }

    long double v = u;
    double x2 = 0.5 * x;
    double f2 = 0.5 * f;
    double f_x_2n = f - x;
    long double t = 0.0, lt = 0.0;

    if (f2 * DBL_EPSILON > 0.125 &&
        fabs(t = x2 - f2) <
        sqrt(DBL_EPSILON) * f2) {

        lt = (1 - t) * (2 - t / (f2 + 1)) - log(sqrt(2 * M_PI)) - 0.5 * log(f2 + 1);
    }
    else {
        /* Usual case 2: careful not to overflow .. : */
        lt = f2 * log(x2) - x2 - SpecialFunctions::Gamma::lgammafn(f2 + 1);
    }

    bool tSml = (lt < _dbl_min_exp);
    l_lam = -1.0;
    double l_x = -1.0;
    double term = 0.0;
    long double ans = 0.0;
    if (tSml) {
        if (x > f + theta + 5 * sqrt(2 * (f + 2 * theta))) {
            return R_DT_1;
        }
        l_x = log(x);
        ans = term = 0.0; t = 0.0;
    }
    else {
        t = expl(lt);
        ans = term = (double)(v * t);
    }

    int temp = 0;
    double f_2n = f + 2.0;
    double bound = 0.0;
    bool is_r = false;
    bool is_b = false;
    f_x_2n += 2.0;
    for (int n = 1; n <= itrmax; n++, f_2n += 2.0, f_x_2n += 2.0) {

//#ifndef MATHLIB_STANDALONE
//        if (n % 1000 == 0) R_CheckUserInterrupt();
//#endif

        if (f_x_2n > 0) {
            bound = (double)(t * x / f_x_2n);
            is_r = false;
            if (((is_b = (bound <= errmax)) &&
                (is_r = (term <= reltol * ans)))) {
                break;
            }
        }

        if (lamSml) {
            lu += l_lam - log(n);
            if (lu >= _dbl_min_exp) {
                v = u = expl(lu);
                lamSml = false;
            }
        }
        else {
            u *= lam / n;
            v += u;
        }
        if (tSml) {
            lt += l_x - log(f_2n);
            if (lt >= _dbl_min_exp) {
                t = expl(lt);
                tSml = false;
            }
        }
        else {
            t *= x / f_2n;
        }
        if (!lamSml && !tSml) {
            term = (double)(v * t);
            ans += term;
        }
        temp = n;
    }

    if (temp > itrmax) {
        std::cout << "Warning cdf(x = " << x << ", f = " << f << ", theta = " << theta
            << ", ...): not converaged in " << itrmax << " iter." << std::endl;
    }

    double dans = (double)ans;
    double tmp1 = log_p ? log1p(-dans) : (0.5 - dans + 0.5);
    double tmp2 = log_p ? log(dans) : dans;
    return lower_tail ? tmp2 : tmp1;
}

double NonCentralChisq::logspaceAdd(double logx, double logy)
{
    return std::max(logx, logy) + log1p(exp(-fabs(logx - logy)));
}

double NonCentralChisq::quantile(double p, double df, double ncp,
    bool lower_tail, bool log_p)
{
    const double accu = 1e-13;
    const double racc = 4 * DBL_EPSILON;
    const double Eps = 1e-11;
    const double rEps = 1e-10;

    double ux = 0.0, lx = 0.0, ux0 = 0.0, nx = 0.0, pp = 0.0;

    if (std::isnan(p) || std::isnan(df) || std::isnan(ncp))
        return p + df + ncp;

    if (!std::isfinite(df)) return InfNaN::nan();
    if (df < 0 || ncp < 0) return InfNaN::nan();

    if (log_p) {
	    if(p > 0)
            return InfNaN::nan();
        if (p == 0)
            return lower_tail ? InfNaN::posinf() : 0.0;
        if (p == InfNaN::neginf())
            return lower_tail ? 0.0 : InfNaN::posinf();
    }
    else {
	    if(p < 0 || p > 1)
            return InfNaN::nan();
	    if(p == 0)
	        return lower_tail ? 0.0 : InfNaN::posinf();
        if (p == 1)
            return lower_tail ? InfNaN::posinf() : 0.0;
    }

    pp = log_p ? exp(p) : p;
    if (pp > 1 - DBL_EPSILON) return lower_tail ? InfNaN::posinf() : 0.0;
    
    {
        double b, c, ff;
        b = (ncp * ncp) / (df + 3 * ncp);
        c = (df + 3 * ncp) / (df + 2 * ncp);
        ff = (df + 2 * ncp) / (c * c);
        ux = b + c * Chisq::quantile(p, ff, lower_tail, log_p);
        if (ux < 0) ux = 1;
        ux0 = ux;
    }

    if (!lower_tail && ncp >= 80) {
        if (pp < 1e-10)
            std::cout << "full precision may not have been achieved in "
                << "Non-central Chisq quantile()" << std::endl;
        p = log_p ? -expm1(p) : (0.5 - (p)+0.5);
        lower_tail = true;
    }
    else {
        p = pp;
    }

    pp = std::min(1 - DBL_EPSILON, p * (1 + Eps));
    if (lower_tail) {
        for (; ux < DBL_MAX &&
            cdfRaw(ux, df, ncp, Eps, rEps, 10000, true, false) < pp;
            ux *= 2);
        pp = p * (1 - Eps);
        for (lx = std::min(ux0, DBL_MAX);
            lx > DBL_MIN &&
            cdfRaw(lx, df, ncp, Eps, rEps, 10000, true, false) > pp;
            lx *= 0.5);
    }
    else {
        for (; ux < DBL_MAX &&
            cdfRaw(ux, df, ncp, Eps, rEps, 10000, false, false) > pp;
            ux *= 2);
        pp = p * (1 - Eps);
        for (lx = std::min(ux0, DBL_MAX);
            lx > DBL_MIN &&
            cdfRaw(lx, df, ncp, Eps, rEps, 10000, false, false) < pp;
            lx *= 0.5);
    }

    if (lower_tail) {
        do {
            nx = 0.5 * (lx + ux);
            if (cdfRaw(nx, df, ncp, accu, racc, 100000, true, false) > p)
                ux = nx;
            else
                lx = nx;
        } while ((ux - lx) / nx > accu);
    }
    else {
        do {
            nx = 0.5 * (lx + ux);
            if (cdfRaw(nx, df, ncp, accu, racc, 100000, false, false) < p)
                ux = nx;
            else
                lx = nx;
        } while ((ux - lx) / nx > accu);
    }
    return 0.5 * (ux + lx);
}

double NonCentralChisq::rand(double df, double lambda)
{
    if (std::isnan(df) || !std::isfinite(lambda) || df < 0.0 || lambda < 0.0) {
        return InfNaN::nan();
    }

    if(lambda == 0.0) {
        return (df == 0.0) ? 0.0 : Gamma::rand(df / 2.0, 2.0);
    }
    else {
        double r = Poisson::rand(lambda / 2.0);
        if (r > 0.0)  r = Chisq::rand(2.0 * r);
        if (df > 0.0) r += Gamma::rand(df / 2.0, 2.0);
        return r;
    }
}
