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

#include "normal.h"

double Normal::dnorm(const double x,
	const double mean, const double sigma)
{
	return gsl_ran_ugaussian_pdf((x - mean) / sigma) / sigma;
}

double Normal::pnorm(const double q,
	const double mean, const double sigma)
{
	return gsl_cdf_ugaussian_P((q - mean) / sigma);
}

double Normal::qnorm(const double p,
	const double mean, const double sigma)
{
	return mean + gsl_cdf_ugaussian_Pinv(p) * sigma;
}

double Normal::rnorm(const double n,
	const double mean, const double sigma)
{
	return 1.0;	// todo
}

/*

class Normal
{
public:
	Normal() = default;
	
	static double dnorm(const double x,
		const double mean = 0.0, const double sigma = 0.0);

	static double pnorm(const double q,
		const double mean = 0.0, const double sigma = 0.0);

	static double qnorm(const double p,
		const double mean = 0.0, const double sigma = 0.0);

	static double rnorm(const double n,
		const double mean = 0.0, const double sigma = 0.0);

	~Normal() = default;
};
*/
