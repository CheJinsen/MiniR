/*
 * Normal distribution functions of density probability and quantile.
 * Copyright (C) 2020 Jinsen Che
 *
 * This file is part of MiniR.
 *
 * MiniR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MiniR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the impiled warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/license/>.
 */

#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>
#include <gsl/gsl_statistics_double.h>


class Base
{
public:
	Base() = default;

	/*
	 * This function returns the arithmetic mean of data.
	 * param data: a vector of double numbers.
	 */
	static double mean(const std::vector<double>& data);

	/*
	 * This function returns the estimated, or sample, variance of data.
	 * param data: a vector of double numbers.
	 */
	static double var(const std::vector<double>& data);

	/*
	 * This functions return the square root of the corresponding
	 * variance function above.
	 */
	static double sd(const std::vector<double>& data);

	/*
	 * This function return the total sum of squares(TSS) fo data about mean.
	 */
	static double tss(const std::vector<double>& data);

	/*
	 * This function computes the skewness of data.
	 */
	static double skew(const std::vector<double>& data);

	/*
	 * This function computes the kurtosis of data.
	 */
	static double kurtosis(const std::vector<double>& data);

	
	~Base() = default;

private:
	static void vector_check(const std::vector<double>& data);
};
