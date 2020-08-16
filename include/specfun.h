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
 * along with MiniR. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include "infnan.h"
#include "acmtoms.h"

namespace SpecialFunctions
{
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