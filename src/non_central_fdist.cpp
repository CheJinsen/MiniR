/*
 * This file is part of MiniR.
 * AUTHOR
 *    Peter Ruckdeschel, peter.ruckdeschel@uni-bayreuth.de.
 *    April 13, 2006.
 *
 * Copyright (C) 1998	Ross Ihaka
 * Copyright (C) 2006-2008 The R Core Team
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

double NonCentralFdist::pdf(double x, double df1, double df2,
    double ncp, bool give_log)
{
    if (std::isnan(x) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
        return x + df2 + df1 + ncp;
    if (df1 <= 0. || df2 <= 0. || ncp < 0)
        return InfNaN::nan();
    if (x < 0.0)
        return give_log ? InfNaN::neginf() : 0.0;
    if (!std::isfinite(ncp))
        return InfNaN::nan();

    if (!std::isfinite(df1) && !std::isfinite(df2)) {
        if (x == 1.0)
            return InfNaN::posinf();
        else
            return give_log ? InfNaN::neginf() : 0.0;
    }
    if (!std::isfinite(df2))
        return df1 * NonCentralChisq::pdf(x * df1, df1, ncp, give_log);
    
    double f = 0.0;
    double z = 0.0;
    if (df1 > 1e14 && ncp < 1e7) {
        f = 1 + ncp / df1;
        z = Gamma::cdf(1.0 / x / f, df2 / 2, 2.0 / df2, give_log);
        return give_log ? z - 2 * log(x) - log(f) : z / (x * x) / f;
    }

    double y = (df1 / df2) * x;
    z = NonCentralBeta::pdf(y / (1 + y), df1 / 2.0, df2 / 2.0, ncp, give_log);
    return  give_log ?
        z + log(df1) - log(df2) - 2 * log1p(y) :
        z * (df1 / df2) / (1 + y) / (1 + y);
}

double NonCentralFdist::cdf(double x, double df1, double df2, double ncp,
    bool lower_tail, bool log_p)
{
    if (std::isnan(x) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
        return x + df2 + df1 + ncp;
    if (df1 <= 0. || df2 <= 0. || ncp < 0)
        return InfNaN::nan();
    if (!std::isfinite(ncp))
        return InfNaN::nan();
    if (!std::isfinite(df1) && !std::isfinite(df2)) /* both +Inf */
        return InfNaN::nan();

    double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
    double R_D__1 = log_p ? 0.0 : 1.0;
    double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
    double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
    if (x <= 0.0) return R_DT_0;
    if (x >= InfNaN::posinf()) return R_DT_1;

    if (df2 > 1e8)
        return NonCentralChisq::cdf(x * df1, df1, ncp, lower_tail, log_p);

    double y = (df1 / df2) * x;
    return Base::nonCentralBetaCdf2(y / (1. + y), 1. / (1. + y), df1 / 2., df2 / 2.,
        ncp, lower_tail, log_p);
}

double NonCentralFdist::quantile(double p, double df1, double df2, double ncp,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(df1) || std::isnan(df2) || std::isnan(ncp))
		return p + df1 + df2 + ncp;
	if (df1 <= 0. || df2 <= 0. || ncp < 0)
		return InfNaN::nan();
	if (!std::isfinite(ncp))
		return InfNaN::nan();
	if (!std::isfinite(df1) && !std::isfinite(df2))
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

	if (df2 > 1e8) /* avoid problems with +Inf and loss of accuracy */
		return NonCentralChisq::quantile(p, df1, ncp, lower_tail, log_p) / df1;

	double y = NonCentralBeta::quantile(p, df1 / 2.0, df2 / 2.0, ncp, lower_tail, log_p);
	return y / (1 - y) * (df2 / df1);
}

double NonCentralFdist::rand(const double df1, const double df2, const double ncp)
{
    return (NonCentralChisq::rand(df1, ncp) / df1) /
		(Chisq::rand(df2) / df2);
}