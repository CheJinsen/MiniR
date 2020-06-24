
#include <iostream>
#include <iomanip>
#include "nmath.h"

int main()
{
	std::cout << std::setprecision(15);
	std::cout << dnorm(0.985) << std::endl;
	std::cout << pnorm(0.975, 0.0, 1.0, true, false) << std::endl;
	std::cout << pnorm(0.975) << std::endl;
	std::cout << pnorm(sqrt(-1))<< std::endl;
	std::cout << pnorm(820, 1019, 209) << std::endl;
	std::cout << (pnorm(1500, 1019, 209) + pnorm(1500, 1019, 209, false))<< std::endl;

	std::cout << qnorm(0.85) << std::endl;

	std::cout << pbinom(2, 6, 0.6) << std::endl;
	std::cout << (1 - pbinom(2, 20, 0.05)) << std::endl;
	std::cout << pbinom(2, 20, 0.05, false) << std::endl;
	std::cout << pbinom(sqrt(-1), 20, 0.05) << std::endl;

	std::cout << dbinom(75, 1500, 0.05) << std::endl;
	std::cout << pbinom(74, 1500, 0.05) << std::endl;

	std::cout << qbinom(0.8, 1500, 0.05) << std::endl;
	std::cout << qbinom(0.9, 1500, 0.05) << std::endl;
	std::cout << qbinom(sqrt(-2), 1500, 0.03) << std::endl;

	std::cout << qnorm(0.85) << std::endl;
	std::cout << qnorm(0.05, 80, 12) << std::endl;
	std::cout << qnorm(0.95, 80, 12) << std::endl;
	std::cout << qnorm(1, 80, 12) << std::endl;

	std::cout << dnorm(1, 2, 3) << std::endl;

	std::cout << pnorm(820, 1019, 209) << std::endl;
	std::cout << pnorm(1500, 1019, 209) << std::endl;
	std::cout << pnorm(2.824) << std::endl;
	std::cout << (pnorm(100, 80, 12) - pnorm(90, 80, 12)) << std::endl;
	std::cout << qnorm(0.975) << std::endl;

	std::cout << qt(0.95, 23) << std::endl;
	std::cout << qt(0.975, 9) << std::endl;

	std::cout << (2.0 * pt(-3, 199)) << std::endl;
	std::cout << "Done..." << std::endl;
	return 0;
}
