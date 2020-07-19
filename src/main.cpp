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

#include <iostream>
#include <iomanip>
#include "randist.h"
#include <random>
#include <ctime>

using namespace Randist;

void describe(const std::vector<double>& v);

// test session
int main(int argc, char* argv[])
{
	std::cout << std::setprecision(15);

	// std::cout << NonCentralTdist::quantile(0.875, 12, 9.87) << std::endl;


#if 1
	const int n = 1000;
	std::vector<double> vec;
	for (int i = 0; i < n; ++i) {
		double rand = NonCentralTdist::rand(12, 9.87);
		std::cout << rand << std::endl;
		vec.push_back(rand);
	}

	std::cout << std::endl;
	describe(vec);
#endif

	return 0;
}

void describe(const std::vector<double>& v)
{
	if (v.empty()) {
		std::cout << "empty vector" << std::endl;
		exit(1);
	}

	double sum = std::accumulate(v.cbegin(), v.cend(), 0.0);
	double mean = sum / v.size();
	double stddev = 0.0;
	double tmp = 0.0;

	for (auto& value : v) {
		tmp += (value - mean) * (value - mean);
	}
	stddev = sqrt(tmp / (v.size() - 1));

	std::cout << "sum       = " << sum << std::endl;
	std::cout << "mean      = " << mean << std::endl;
	std::cout << "stddev    = " << stddev << std::endl;
}
