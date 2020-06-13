/*
 * Binomial distribution functions of density probability and quantile.
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


class Binomial
{
public:
	Binomial() = default;
	
	/*
	 * This function compute the probability p(k) of obtaining
	 * k from a binomial distribution with parameters size and prob.
	 * param x: quantiles of double.
	 * param size: number of trials (zero or more).
	 * param prob: probability of success on each trial.
	 */
	static double dbinom(const double x,
		const double size, const double prob);
	
	/*
	 * This function compute the cumulative distribution functions
	 * for the binomial distribution with parameters size and prob.
	 */
	static double pbinom(const double q,
		const double size, const double prob);

	/*
	 * This function compute the inverse of cumulative distribution.
	 */
	static double qbinom(const double p,
		const double size, const double prob);

	/*
	 * This function returns a rnadom integer from the binomial
	 * distribution, the number of successes in n independent trials
	 * with probability p.
	 */
	static std::vector<double> rbinom(const double n,
		const double size, const double prob);
	

	~Binomial() = default;

private:
	static double do_search(double y, double *z, double p,
		double n, double pr, double incr);
};
