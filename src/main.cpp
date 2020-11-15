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

#include <vector>
#include <iostream>
#include <iomanip>
#include "regression.h"
#include <Eigen/Core>

using namespace Regression;

int main(int argc, char* argv[])
{
	std::cout << std::setprecision(15);

	const std::vector<double> estriol = {7,9,9,12,14,16,16,14,16,16,17,19,21,24,15,16,17,25,27,15,
		15,15,16,19,18,17,18,20,22,25,24};
	const std::vector<double> birthweight = {25,25,25,27,27,27,24,30,30,31,30,31,30,28,32,32,32,32,34,
		34,34,35,35,34,35,36,37,38,40,39,43};

	Regs::lm(birthweight, estriol);

	Eigen::MatrixXcf a = Eigen::MatrixXcf::Random(3, 3);
	std::cout << a << std::endl;

	Eigen::MatrixXcf m(3, 3);
	m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	std::cout << "m = \n";
	std::cout << m << "\n";
	return 0;
}
