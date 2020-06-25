
#include <iostream>
#include <iomanip>
#include "nmath.h"

int main()
{
	std::cout << std::setprecision(15);
	std::cout << dbeta(0.2, 0.3, 0.4) << std::endl;
	std::cout << qbeta(0.2, 0.3, 0.4) << std::endl;
	return 0;
}

