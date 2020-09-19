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

#include "statistics.h"

using namespace Statistics;

double Stats::mean(const std::vector<double>& v)
{
	if (v.empty()) {
		std::cout << "vector is empty." << std::endl;
		return 0.0;
	}
	
	double sum = accumulate(v.cbegin(), v.cend(), 0.0);
	return sum / v.size();
}