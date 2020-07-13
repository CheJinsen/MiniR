#include <iostream>
#include <iomanip>
#include "randist.h"

using namespace Randist;

// test session
int main(int argc, char* argv[])
{
	std::cout << std::setprecision(15);

	std::cout << Tukey::cdf(1, 1, 2, 3) << std::endl;
	std::cout << Tukey::quantile(0.95, 1, 2, 3) << std::endl;
	std::cout << Wilcox::pdf(23, 20, 3) << std::endl;
	std::cout << Wilcox::cdf(151, 205, 3) << std::endl;
	std::cout << Wilcox::quantile(0.912, 13, 200) << std::endl;

	std::cout << Signrank::pdf(100, 200, true) << std::endl;
	std::cout << Signrank::cdf(0.975, 10) << std::endl;
	std::cout << Signrank::quantile(0.975, 190) << std::endl;
	return 0;
}
