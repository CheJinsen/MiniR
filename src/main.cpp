
/*
 * test session
 */

#include <iostream>
#include <iomanip>
#include "normal.h"

int main(int argc, char *argv[])
{
	std::cout << std::setprecision(18);
	std::cout << Normal::dnorm(2, 3, 4) << std::endl;
	std::cout << Normal::pnorm(0.975) << std::endl;
	std::cout << (Normal::pnorm(100, 80, 12) - Normal::pnorm(90, 80, 12));
	std::cout << std::endl;

	std::cout << Normal::qnorm(0.85) << std::endl;
	std::cout << Normal::qnorm(0.05, 80, 12) << std::endl;
	std::cout << Normal::qnorm(0.05) << std::endl;
	std::cout << Normal::qnorm(0.95, 80, 12) << std::endl;
	std::cout << "Done..." << std::endl;
}
