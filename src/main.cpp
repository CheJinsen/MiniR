#include <iostream>
#include <iomanip>
#include "randist.h"
#include <random>
#include <ctime>

using namespace Randist;

void describe(const std::vector<double>& v);

// test session
int main(int argc, char* argv[])
{
	std::cout << std::setprecision(15);

	const int n = 1000;
	std::vector<double> vec;
	for (int i = 0; i < n; ++i) {
		double rand = Beta::rand();
		std::cout << rand << std::endl;
		vec.push_back(rand);
	}

	std::cout << std::endl;
	describe(vec);

	return 0;
}

void describe(const std::vector<double>& v)
{
	if (v.empty()) {
		std::cout << "empty vector" << std::endl;
		exit(1);
	}

	double sum = std::accumulate(v.cbegin(), v.cend(), 0.0);
	double mean = sum / v.size();
	double stddev = 0.0;
	double tmp = 0.0;

	for (auto& value : v) {
		tmp += (value - mean) * (value - mean);
	}
	stddev = sqrt(tmp / (v.size() - 1));

	std::cout << "sum       = " << sum << std::endl;
	std::cout << "mean      = " << mean << std::endl;
	std::cout << "stddev    = " << stddev << std::endl;
}
