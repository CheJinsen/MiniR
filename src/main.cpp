
#include <iostream>
#include <iomanip>
#include "randist.h"

int main()
{
	std::cout << std::setprecision(15);
	std::cout << dnorm(0.985) << std::endl;
	std::cout << "Done..." << std::endl;
	return 0;
}
