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

#include "base.h"

double Base::mean(const std::vector<double>& data)
{
	vector_check(data);
	double d[data.size()];
	std::copy(data.cbegin(), data.cend(), d);

	return gsl_stats_mean(d, 1, data.size());
}

double Base::var(const std::vector<double>& data)
{
	vector_check(data);
	double d[data.size()];
	std::copy(data.cbegin(), data.cend(), d);

	return gsl_stats_variance(d, 1, data.size());
}

double Base::sd(const std::vector<double>& data)
{
	vector_check(data);
	double d[data.size()];
	std::copy(data.cbegin(), data.cend(), d);

	return gsl_stats_sd(d, 1, data.size());
}

double Base::tss(const std::vector<double>& data)
{
	vector_check(data);
	double d[data.size()];
	std::copy(data.cbegin(), data.cend(), d);

	return gsl_stats_tss(d, 1, data.size());
}

void Base::vector_check(const std::vector<double>& data)
{
	if (data.empty()) {
		std::cerr << "Vector is empty." << std::endl;
		exit(1);
	}
}

