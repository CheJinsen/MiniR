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

#include "infnan.h"

double InfNaN::divf(const double a, const double b)
{
	return a / b;
}

double InfNaN::posinf()
{
	return divf(1.0, 0.0);
}

double InfNaN::neginf()
{
	return divf(-1.0, 0.0);
}

double InfNaN::nan()
{
	return divf(0.0, 0.0);
}