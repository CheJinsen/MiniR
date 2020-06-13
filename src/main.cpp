
/*
 * test session
 */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <numeric>
#include "normal.h"
#include "base.h"

int main(int argc, char *argv[])
{
	std::cout << std::setprecision(10);
	std::vector<double> vecd = {1.1, 2.1, 3.2, 2.2, 5.01, 4.4};
	std::cout << "mean = " << Base::mean(vecd) << std::endl;
	std::cout << "variance = " << Base::var(vecd) << std::endl;
	std::cout << "sd = " << Base::sd(vecd) << std::endl;
	std::cout << "tss = " << Base::tss(vecd) << std::endl;
	return 0;
}




int main01(int argc, char *argv[])
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

	return 0;
}
