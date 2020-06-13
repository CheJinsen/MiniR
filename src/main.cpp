
/*
 * test session
 */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <numeric>
#include "normal.h"

int main(int argc, char *argv[])
{
	std::cout << std::setprecision(15);
	std::cout << Normal::dnorm(2, 3, 4) << std::endl;
	std::cout << Normal::pnorm(0.975) << std::endl;
	std::cout << (Normal::pnorm(100, 80, 12) - Normal::pnorm(90, 80, 12));
	std::cout << std::endl;

	std::cout << Normal::qnorm(0.85) << std::endl;
	std::cout << Normal::qnorm(0.05, 80, 12) << std::endl;
	std::cout << Normal::qnorm(0.05) << std::endl;
	std::cout << Normal::qnorm(0.95, 80, 12) << std::endl;
	
	std::vector<double> ran;
	ran = Normal::rnorm(1000);
	double sum = std::accumulate(ran.cbegin(), ran.cend(), 0.0);
	std::cout << "sum = " << sum << std::endl;
	std::cout << "mean = " << sum / ran.size() << std::endl;

	Normal::set_seed(1234567890);
	ran = Normal::rnorm(1000);
	std::cout << "now sum = "
		<< std::accumulate(ran.cbegin(), ran.cend(), 0.0) << std::endl;
	std::cout << "Done..." << std::endl;
}
