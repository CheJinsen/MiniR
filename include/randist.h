/*
 * This file is part of MiniR.
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
 * along with Foobar. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cfloat>
#include <vector>
#include <random>
#include <ctime>
#include <cstdlib>

#include "basedist.h"
#include "infnan.h"
#include "acmtoms.h"
#include "specfun.h"

namespace Randist
{
	class Normal	// Normal Distribution
	{
	public:
		static double pdf(double x, double mean = 0.0,
			double sd = 1.0, bool log_p = false);
		static double cdf(double q, double mean = 0.0,
			double sd = 1.0, bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double mean = 0.0,
			double sd = 1.0, bool lower_tail = true, bool log_p = false);
		static double rand(const double mean = 0.0, const double sd = 1.0);

	private:
		static void cdfBoth(double x, double& cum, double& ccum, bool i_tail, bool log_p);
	};

	class Uniform 	// Uniform Distribution
	{
	public:
		static double pdf(double x, double min = 0.0,
			double max = 1, bool give_log = false);
		static double cdf(double q, double min = 0.0,
			double max = 1, bool lower_tail = true, bool give_log = false);
		static double quantile(double x, double min = 0.0,
			double max = 1, bool lower_tail = true, bool give_log = false);
		static double rand(const double min = 0.0, const double max = 1.0);
	};

	class Gamma : private Base 	// Gamma Distributon
	{
	public:
		static double pdf(double x, double shape,
			double scale, bool give_log = false);
		static double cdf(double q, double shape,
			double scale, bool lower_tail = true, bool give_log = false);
		static double quantile(double q, double shape,
			double scale, bool lower_tail = true, bool give_log = false);
		static double rand(const double shape = 1.0, const double scale = 1.0);

	private:
		static double logcf(double x, double i, double d, double eps);
		static inline double sqr(double x) { return x * x; }
		static double log1pmx(double x);
		static double lgamma1p(double a);
		static double logspaceAdd(double logx, double logy);
		static double logspaceSub(double logx, double logy);
		static double logspaceSum(const double* logx, int n);
		static double poissonPdfWrap(double x_plus_1, double lambda, bool give_log);
		static double cdfSmallx(double x, double alph, bool lower_tail, bool log_p);
		static double pdUpperSeries(double x, double y, bool log_p);
		static double pdLowercf(double y, double d);
		static double pdLowerSeries(double lambda, double y);
		static double dpnorm(double x, bool lower_tail, double lp);
		static double poissonCdfAsymp(double x, double lambda, bool lower_tail, bool log_p);
		static double cdfRaw(double x, double alph, bool lower_tail, bool log_p);
		static double chisqQuantileAppr(double p, double nu, double g,
			bool lower_tail, bool log_p, double tol);
	};

	class Beta : private Base 	// Beta Distribution
	{
	public:
		static double pdf(double x, double a, double b, bool log_p = false);
		static double cdf(double q, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double shape1, const double shape2);

	private:
		static double cdfRaw(double x, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static void quantileRaw(double alpha, double p, double q,
			bool lower_tail, bool log_p, int swap_01, double log_q_cut, int n_N, double* qb);
		static double powDi(double x, int i);
		static void showBet(double aa, double beta, double u1, double& v, double& w);
	};

	class Lognormal 	// Lognormal Distribution
	{
	public:
		static double pdf(double x, double meanlog, double sdlog,
			bool give_log = false);
		static double cdf(double q, double meanlog, double sdlog,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double meanlog, double sdlog,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double meanlog = 0.0, const double sdlog = 1.0);
	};

	class Chisq 	// Chi-squared Distribution
	{
	public:
		static double pdf(double x, double df, bool give_log = false);
		static double cdf(double q, double df,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double df,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double df = 1.0);
	};

	class NonCentralChisq : private Base 	// Non-central Chi-squared Distribution
	{
	public:
		static double pdf(double x, double df, double ncp,
			bool give_log = false);
		static double cdf(double q, double df, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double df, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double df, const double lambda);

	private:
		static double cdfRaw(double x, double f, double theta,
			double errmax, double reltol, int itrmax,
			bool lower_tail, bool log_p);
		static double logspaceAdd(double logx, double logy);
	};

	class Fdist :private Base	// F Distribution
	{
	public:
		static double pdf(double x, double df1, double df2,
			bool give_log = false);
		static double cdf(double x, double df1, double df2,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double x, double df1, double df2,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double df1, const double df2);
	};

	class Tdist : private Base	// Student t Distribution
	{
	public:
		static double pdf(double x, double df, bool log_p = false);
		static double cdf(double q, double df,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double df,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double df);
	};

	class Binomial : private Base 	// Binomial Distribution
	{
	public:
		static double pdf(double x, double size, double prob, bool log_p = false);
		static double cdf(double q, double size, double prob,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double size, double prob,
			bool lower_tail = true, bool log_p = false);
		static int rand(const int size, const double prob);

	private:
		static double doSearch(double y, double* z, double p,
			double n, double pr, double incr);
	};

	class Cauchy : private Base 	// Cauchy Distribution
	{
	public:
		static double pdf(double x, double location = 0.0,
			double scale = 1.0, bool give_log = false);
		static double cdf(double q, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double location = 0.0, const double scale = 1.0);
	};

	class Exp	// Exponential Distribution
	{
	public:
		static double pdf(double x, double rate = 1.0, bool log_p = false);
		static double cdf(double q, double rate = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double rate = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double lambda = 1.0);
	};

	class Geom : private Base	// Geometric Distribution
	{
	public:
		static double pdf(double x, double prob, bool log_p = false);
		static double cdf(double q, double prob,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double prob,
			bool lower_tail = true, bool log_p = false);
		static int rand(const double prob = 0.5);
	};

	class Hyper : private Base	// Hypergeometric Distribution
	{
	public:
		static double pdf(double x, double m, double n,
			double k, bool log_p = false);
		static double cdf(double q, double m, double n,
			double k, bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double m, double n,
			double k, bool lower_tail = true, bool log_p = false);
		static double rand(double m, double n, double k);

	private:
		static bool isNegOrNonInt(double x);	// negative or non integer
		static double pdhyper(double x, double NR, double NB,
			double n, bool log_p);
		static double afc(int i);
	};

	class NegBinomial : private Base	// Negative Binomial Distribution
	{
	public:
		static double pdf(double x, double size, double prob, bool log_p = false);
		static double cdf(double q, double size, double prob,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double size, double prob,
			bool lower_tail = true, bool log_p = false);

	
		static double pdf_mu(double x, double size, double mu, bool log_p = false);
		static double cdf_mu(double q, double size, double mu,
			bool lower_tail = true, bool log_p = false);
		static double quantile_mu(double p, double size, double mu,
			bool lower_tail = true, bool log_p = false);
		static int rand(const int size = 1, const double prob = 0.5);

	private:
		static double doSearch(double y, double* z, double p,
			double n, double pr, double incr);
	};

	class Poisson : private Base 	// Poisson Distribution
	{
	public:
		static double pdf(double x, double lambda, bool give_log = false);
		static double cdf(double q, double lambda,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double lambda,
			bool lower_tail = true, bool give_log = false);
		static int rand(const double lambda);

	private:
		static double doSearch(double y, double* z, double p,
			double lambda, double incr);
	};

	class Weibull 	// Weibull Distribution
	{
	public:
		static double pdf(double x, double shape, double scale = 1.0,
			bool give_log = false);
		static double cdf(double q, double shape, double scale = 1.0,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double shape, double scale = 1.0,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double shape = 1.0, const double scale = 1.0);
	};

	class Logistic 	// Logistic Distribution
	{
	public:
		static double pdf(double x, double location = 0.0, double scale = 1.0,
			bool log_p = false);
		static double cdf(double q, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double location = 0.0, const double scale = 1.0);

	private:
		static double log1pexp(const double x);
	};

	class NonCentralBeta : private Base 	// Non-central Beta Distribution
	{
	public:
		static double pdf(double x, double a, double b, double ncp,
			bool give_log = false);
		static double cdf(double q, double a, double b, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double a, double b, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double rand(const double shape1, const double shape2, const double ncp = 0.0);
	};

	class NonCentralFdist : private Base	// Non-central F Distribution
	{
	public:
		static double pdf(double x, double df1, double df2,
			double ncp, bool give_log = false);
		static double cdf(double x, double df1, double df2,
			double ncp, bool lower_tail = true, bool give_log = false);
		static double quantile(double x, double df1, double df2,
			double ncp, bool lower_tail = true, bool give_log = false);
		static double rand(const double df1, const double df2, const double ncp);
	};

	class NonCentralTdist	// Non-Central Student t Distribution
	{
	public:
		static double pdf(double x, double df, double ncp, bool log_p = false);
		static double cdf(double q, double df, double ncp,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double df, double ncp,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double df, const double ncp);
	};

	class Tukey	// Studentized Range Distribution
	{
	public:
		static double cdf(double q, double nranges, double nmeans,
			double df, bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double nranges, double nmeans,
			double df, bool lower_tail = true, bool log_p = false);

	private:
		static double wprob(double w, double rr, double cc);
		static double qinv(double p, double c, double v);
		static double dtvalue(const double x,
			const bool lower_tail, const bool log_p);
	};

	class Wilcox	// Wilcoxon Rank Sum Distribution
	{
	public:
		static double pdf(double x, double m, double n, bool log_p = false);
		static double cdf(double q, double m, double n,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double m, double n,
			bool lower_tail = true, bool log_p = false);
		static double rand(double m, double n);

	private:
		static int allocated_m;
		static int allocated_n;
		static std::vector<std::vector<std::vector<double>>> w;

		static void wInitMaybe(int m, int n);
		static double cwilcox(int k, int m, int n);
		static double dtvalue(const double x, const bool lower_tail, const bool log_p);
		static double rbits(const int bits);
		static double uniformIndex(const double dn);
	};

	class Signrank	// Wilcoxon Signed Rank Distribution
	{
	public:
		static double pdf(double x, double n, bool log_p = false);
		static double cdf(double q, double n,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double n,
			bool lower_tail = true, bool log_p = false);
		static double rand(double n);

	private:
		static std::vector<double> w;
		static int allocated_n;

		static void wInitMaybe(int n);
		static double csignrank(int k, int n);
		static double dtvalue(const double x,
			const bool lower_tail, const bool log_p);
	};
}