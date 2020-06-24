/*
 * sys/infnan.c
 *
 * Copyright (C) 2001, 2004, 2007, 2010 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "nmath.h"

double fdiv(const double a, const double b)
{
	return a / b;
}

double posinf()
{
	return fdiv(1.0, 0.0);
}

double neginf()
{
	return fdiv(-1.0, 0.0);
}

double nan()
{
	return fdiv(0.0, 0.0);
}