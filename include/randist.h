#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cfloat>
#include <vector>
#include <random>
#include <ctime>
#include <cstdlib>

#include "infnan.h"
#include "acmtoms.h"
#include "specfun.h"

namespace Randist
{
	constexpr double M_LN_SQRT_2PI	= 0.918938533204672741780329736406;	// log(sqrt(2*pi))
	constexpr double M_LN_2PI		= 1.837877066409345483560659472811;	// log(2*pi)
	// constexpr double M_PI			= 3.141592653589793238462643383280;
	constexpr double M_1_SQRT_2PI	= 0.398942280401432677939946059934;	// 1/sqrt(2pi)
	// constexpr double M_LN2			= 0.693147180559945309417232121458;	// ln(2)
	constexpr double M_SQRT_32		= 5.656854249492380195206754896838;	// sqrt(32)
	constexpr double M_2PI			= 6.283185307179586476925286766559;
	// constexpr double M_LN10			= 2.302585092994045684017991454684; // ln(10)
	// constexpr double M_1_PI			= 0.318309886183790671537767526745;	// 1/pi
	// constexpr double M_PI_2			= 1.570796326794896619231321691640; // pi/2


	class Normal
	{
	public:
		static double pdf(double x, double mean = 0.0,
			double sd = 1.0, bool log_p = false);
		static double cdf(double q, double mean = 0.0,
			double sd = 1.0, bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double mean = 0.0,
			double sd = 1.0, bool lower_tail = true, bool log_p = false);
		static double rand(const double mean = 0.0, const double sd = 1.0);
		static void setSeed(const unsigned long int s);

	private:
		static unsigned long int seed;
		static void cdfBoth(double x, double& cum, double& ccum, bool i_tail, bool log_p);
	};

	class Uniform
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

	class Gamma
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
		static double poissonPdfRaw(double x, double lambda, bool give_log);
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);

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

	class Beta
	{
	public:
		static double pdf(double x, double a, double b, bool log_p = false);
		static double cdf(double q, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static double rand(const double shape1 = 1.0, const double scale = 1.0);
		// static double rand(const double shape1 = 1.0,
		// 	const double scale = 1.0, const double ncp = 0.0);	// todo

	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double cdfRaw(double x, double a, double b,
			bool lower_tail = true, bool log_p = false);
		static double binomialPdfRaw(const double x, const double n,
			const double p, const double q, bool give_lop);
		static void quantileRaw(double alpha, double p, double q,
			bool lower_tail, bool log_p, int swap_01, double log_q_cut,
			int n_N, double* qb);
		static double powDi(double x, int i);
		static void showBet(double aa, double beta, double u1, double& v, double& w);
	};

	class Lognormal
	{
	public:
		static double pdf(double x, double meanlog, double sdlog,
			bool give_log = false);
		static double cdf(double q, double meanlog, double sdlog,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double meanlog, double sdlog,
			bool lower_tail = true, bool give_log = false);
	};

	class Chisq
	{
	public:
		static double pdf(double x, double df, bool give_log = false);
		static double cdf(double q, double df,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double df,
			bool lower_tail = true, bool give_log = false);
	};

	class NonCentralChisq
	{
	public:
		static double pdf(double x, double df, double ncp,
			bool give_log = false);
		static double cdf(double q, double df, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double df, double ncp,
			bool lower_tail = true, bool give_log = false);

	private:
		static double poissonPdfRaw(double x, double lambda, bool give_log);
		static double bd0(const double x, const double np);
		static double stirlerr(const double n);
		static double cdfRaw(double x, double f, double theta,
			double errmax, double reltol, int itrmax,
			bool lower_tail, bool log_p);
		static double logspaceAdd(double logx, double logy);
	};

	class Fdist	// F Distribution
	{
	public:
		static double pdf(double x, double df1, double df2,
			bool give_log = false);
		static double cdf(double x, double df1, double df2,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double x, double df1, double df2,
			bool lower_tail = true, bool give_log = false);
	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double BinomialPdfRaw(const double x, const double n,
			const double p, const double q, bool give_lop);
	};

	class Tdist	// Student t Distribution
	{
	public:
		static double pdf(double x, double df, bool log_p = false);
		static double cdf(double q, double df,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double df,
			bool lower_tail = true, bool log_p = false);
		//double rt(double);
		
	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double tanpi(double x);
	};

	class Binomial
	{
	public:

		static double pdf(double x, double size, double prob, bool log_p = false);
		static double cdf(double q, double size, double prob,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double size, double prob,
			bool lower_tail = true, bool log_p = false);

	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double pdfRaw(const double x, const double n, const double p,
			const double q, bool give_lop);
		static double doSearch(double y, double* z, double p,
			double n, double pr, double incr);
	};

	class Cauchy
	{
	public:
		static double pdf(double x, double location = 0.0,
			double scale = 1.0, bool give_log = false);
		static double cdf(double q, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);

	private:
		static double tanpi(double x);
	};

	class Exp	// Exponential Distribution
	{
	public:
		static double pdf(double x, double rate = 1.0, bool log_p = false);
		static double cdf(double q, double rate = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double rate = 1.0,
			bool lower_tail = true, bool log_p = false);
	};

	class Geom	// Geometric Distribution
	{
	public:
		static double pdf(double x, double prob, bool log_p = false);
		static double cdf(double q, double prob,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double prob,
			bool lower_tail = true, bool log_p = false);

	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double BinomialPdfRaw(const double x, const double n,
			const double p, const double q, bool give_lop);
	};

	class Hyper	// Hypergeometric Distribution
	{
	public:
		static double pdf(double x, double m, double n,
			double k, bool log_p = false);
		static double cdf(double q, double m, double n,
			double k, bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double m, double n,
			double k, bool lower_tail = true, bool log_p = false);

	private:
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double BinomialPdfRaw(const double x, const double n,
			const double p, const double q, bool give_lop);
		static bool isNegOrNonInt(double x);	// negative or non integer
		static double pdhyper(double x, double NR, double NB,
			double n, bool log_p);
	};

	class NegBinomial	// Negative Binomial Distribution
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

	private:
		static double poissonPdfRaw(double x, double lambda, bool give_log);
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double BinomialPdfRaw(const double x, const double n,
			const double p, const double q, bool give_lop);
		static double doSearch(double y, double* z, double p,
			double n, double pr, double incr);
	};

	class Poisson
	{
	public:
		static double pdf(double x, double lambda, bool give_log = false);
		static double cdf(double q, double lambda,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double lambda,
			bool lower_tail = true, bool give_log = false);

	private:
		static double pdfRaw(double x, double lambda, bool give_log);
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static double doSearch(double y, double* z, double p,
			double lambda, double incr);
	};

	class Weibull
	{
	public:
		static double pdf(double x, double shape, double scale = 1.0,
			bool give_log = false);
		static double cdf(double q, double shape, double scale = 1.0,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double shape, double scale = 1.0,
			bool lower_tail = true, bool give_log = false);
	};

	class Logistic
	{
	public:
		static double pdf(double x, double location = 0.0, double scale = 1.0,
			bool log_p = false);
		static double cdf(double q, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double location = 0.0, double scale = 1.0,
			bool lower_tail = true, bool log_p = false);

	private:
		static double log1pexp(const double x);
	};

	class NonCentralBeta
	{
	public:
		static double pdf(double x, double a, double b, double ncp,
			bool give_log = false);
		static double cdf(double q, double a, double b, double ncp,
			bool lower_tail = true, bool give_log = false);
		static double quantile(double p, double a, double b, double ncp,
			bool lower_tail = true, bool give_log = false);

	private:
		static double poissonPdfRaw(double x, double lambda, bool give_log);
		static double bd0(const double x, const double np);
		static double stirlerr(const double x);
		static long double cdfRaw(double x, double o_x, double a,
			double b, double ncp);
		static double cdf2(double x, double o_x, double a, double b, double ncp,
			bool lower_tail, bool log_p);
	};

	class NonCentralFdist	// F Distribution
	{
	public:
		static double pdf(double x, double df1, double df2,
			double ncp, bool give_log = false);
		static double cdf(double x, double df1, double df2,
			double ncp, bool lower_tail = true, bool give_log = false);
		static double quantile(double x, double df1, double df2,
			double ncp, bool lower_tail = true, bool give_log = false);

	private:
		static long double NonCentralBetaCdfRaw(double x, double o_x, double a,
			double b, double ncp);
		static double NonCentralBetaCdf2(double x, double o_x, double a, double b, double ncp,
			bool lower_tail, bool log_p);
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

		// may be have bugs, when n is two big, it's two slow
		// example: quantile(0.912, 13, 2000); ==> 15840
		static double quantile(double p, double m, double n,
			bool lower_tail = true, bool log_p = false);

	private:
		static int allocated_m;
		static int allocated_n;
		static std::vector<std::vector<std::vector<double>>> w;

		static void wInitMaybe(int m, int n);
		static double cwilcox(int k, int m, int n);
		static double dtvalue(const double x,
			const bool lower_tail, const bool log_p);
	};

	class Signrank	// Wilcoxon Signed Rank Distribution
	{
	public:
		static double pdf(double x, double n, bool log_p = false);
		static double cdf(double q, double n,
			bool lower_tail = true, bool log_p = false);
		static double quantile(double p, double n,
			bool lower_tail = true, bool log_p = false);
		// double rsignrank(double);

	private:
		static std::vector<double> w;
		static int allocated_n;

		static void wInitMaybe(int n);
		static double csignrank(int k, int n);
		static double dtvalue(const double x,
			const bool lower_tail, const bool log_p);
	};
}