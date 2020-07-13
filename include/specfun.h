#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include "infnan.h"
#include "acmtoms.h"

namespace SpecialFunctions
{
	const double M_LN_SQRT_PId2 = 0.225791352644727432363097614947;	// log(sqrt(pi/2))
	const double M_LN_SQRT_2PI = 0.918938533204672741780329736406;	// log(sqrt(2*pi))
	const double M_LN_2PI = 1.837877066409345483560659472811;	// log(2*pi)
	// const double M_PI = 3.141592653589793238462643383280;
	const double M_1_SQRT_2PI = 0.398942280401432677939946059934;	// 1/sqrt(2pi)
	// const double M_LN2 = 0.693147180559945309417232121458;	/* ln(2) */
	const double M_SQRT_32 = 5.656854249492380195206754896838;	// sqrt(32)
	const double M_SQRT_PI = 1.772453850905516027298167483341;	// sqrt(pi)

	class Gamma
	{
	public:
		static double gammafn(const double x);
		static double lgammafn(const double x);
		static double lgammafnSign(const double x, int* sgn);
		static double lgammacor(const double x);

		static double logspaceAdd(const double logx, const double logy);

	private:
		static double stirlerr(const double x);
		static int chebyshevInit(const double* dos,
			const int nos, const double eta);
		static double chebyshevEval(const double x,
			const double* a, const int n);
		static double sinpi(const double x);
	};

	class Beta
	{
	public:
		static double lbeta(const double a, const double b);

	};

	class Choose
	{
	public:
		static double choose(double n, double k);
		static double lchoose(double n, double k);
		static double lfastchoose(double n, double k);
		static double lfastchoose2(double n, double k, int* s_choose);

	private:
		static bool isOdd(const double k);
		static bool isInt(const double x);
	};
}