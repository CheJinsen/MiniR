#pragma once

class InfNaN
{
public:
	static double posinf();	// 1.0 / 0.0
	static double neginf();	// -1.0 / 0.0
	static double nan();
	static double divf(const double a, const double b);
};