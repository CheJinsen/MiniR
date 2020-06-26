
#include <iostream>
#include <iomanip>
#include "nmath.h"

int main(int argc, char *argv[])
{
	std::cout << std::setprecision(15);
	std::cout << dunif(0.3, 0.2, 0.5) << std::endl;
	std::cout << dlnorm(3.24) << std::endl;
	std::cout << plnorm(7.23) << std::endl;
	std::cout << qlnorm(0.975) << std::endl;

	std::cout << qf(0.05, 6, 8) << std::endl;
	std::cout << qf(0.99, 5, 9) << std::endl;

	std::cout << dcauchy(1, 2, 3) << std::endl;
	std::cout << pcauchy(3, 2, 1) << std::endl;
	std::cout << qcauchy(0.43, 3, 4) << std::endl;
	return 0;
}
