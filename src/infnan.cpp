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