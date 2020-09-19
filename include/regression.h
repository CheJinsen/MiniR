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
#include <vector>
#include <cstdlib>
#include <iterator>
#include <iostream>
#include "statistics.h"

namespace Regression
{
	class Regs
	{
	public:
		static void lm(const std::vector<double>& y, const std::vector<double>& x);

	private:
		static double correctedSumSquares(const std::vector<double>& v);
		static double rawSumCrossProducts(const std::vector<double>& v1,
			const std::vector<double>& v2);
	};
}