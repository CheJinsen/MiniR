
#include <iostream>
#include <iomanip>
#include "nmath.h"

int main(int argc, char *argv[])
{
	std::cout << choose(20000, 10) << std::endl;
	std::cout << choose(20000, 100) << std::endl;
	std::cout << choose(20000, 1000) << std::endl;

	std::cout << dsignrank(100, 20) << std::endl;
	std::cout << psignrank(9, 5) << std::endl;
	return 0;
}
