
#include <iostream>
#include <iomanip>
#include "nmath.h"

int main()
{
	std::cout << std::setprecision(15);
	std::cout << dbeta(0.2, 0.3, 0.4) << std::endl;
	std::cout << qbeta(0.2, 0.3, 0.4) << std::endl;
	std::cout << dpois(1, 2.3) << std::endl;
	std::cout << dchisq(15, 10) << std::endl;
	std::cout << pchisq(20.48, 10) << std::endl;

	std::cout << qchisq(0.975, 10) << std::endl;
	std::cout << qchisq(0.025, 10) << std::endl;
	std::cout << (1 - ppois(3, 1.15)) << std::endl;
	std::cout << qpois(0.97, 2.3) << std::endl;
	return 0;
}

