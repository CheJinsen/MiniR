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

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


class Normal
{
public:
	Normal() = default;
	
	/*
	 * This function compute the probability density p(x) at x for
	 * a Gaussian distribution with standard deviation mean = 0.0
	 * and sigma = 1.0 default.
	 */
	static double dnorm(const double x,
		const double mean = 0.0, const double sigma = 1.0);
	
	/*
	 * This function compute the cumulative distribution functions.
	 * mean = 0.0 and sigma = 1.0 default.
	 */
	static double pnorm(const double q,
		const double mean = 0.0, const double sigma = 1.0);
	
	/*
	 * This function compute the inverses cumulative distribution
	 * function p(x), mean = 0.0, sigma = 1.0 default.
	 */
	static double qnorm(const double p,
		const double mean = 0.0, const double sigma = 1.0);
	
	/*
	 * This function returns a Gaussian random variate,
	 * with mean zero and standard deviation sigma.
	 * param n: number of observations.
	 * param mean: means of random number.
	 * param sigma: sigma of random number.
	 * return a vector of double random number.
	 */
	static std::vector<double> rnorm(const double n,
		const double mean = 0.0, const double sigma = 1.0);

	static void set_seed(unsigned long int seed);

	~Normal() = default;

private:
	static unsigned long int seed; 
};
