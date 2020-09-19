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

#include "regression.h"

using namespace Regression;
using namespace Statistics;

double Regs::correctedSumSquares(const std::vector<double>& v)
{
	if (v.empty()) {
		std::cout << "vector is empty." << std::endl;
		return 0.0;
	}

	double sum = 0.0;
	double m = Stats::mean(v);
	for (auto& i : v) {
		sum += (i - m) * (i - m);
	}
	return sum;
}

double Regs::rawSumCrossProducts(const std::vector<double>& v1,
	const std::vector<double>& v2)
{
	if (v1.empty() || v2.empty()) {
		std::cout << "vector is empty." << std::endl;
		return 0.0;
	}

	if (v1.size() != v2.size()) {
		std::cout << "difference length of two vectors" << std::endl;
		exit(0);
	}

	double m_v1 = Stats::mean(v1);
	double m_v2 = Stats::mean(v2);
	double sum = 0.0;

	std::vector<double>::const_iterator it_1 = v1.cbegin();
	std::vector<double>::const_iterator it_2 = v2.cbegin();
	while (it_1 != v1.cend() && it_2 != v2.cend()) {
		sum += (*it_1 - m_v1) * (*it_2 - m_v2);
		++it_1;
		++it_2;
	}

	return sum;
}

void Regs::lm(const std::vector<double>& y, const std::vector<double>& x)
{
	double lyx = rawSumCrossProducts(y, x);
	double lxx = correctedSumSquares(x);
	double b = lyx / lxx;
	double intercept = Stats::mean(y) - b * Stats::mean(x);

	std::cout << "intercept = " << intercept << std::endl;
	std::cout << "slop = " << b << std::endl;
}