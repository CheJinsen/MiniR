
#include <iostream>
#include <iomanip>
#include "randist.h"

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
	std::cout << "Done..." << std::endl;
	return 0;
}
