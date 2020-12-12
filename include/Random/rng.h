
#pragma once

#include <random>

namespace minilab
{
namespace random
{

class RNG
{
public:
	RNG();
	RNG(unsigned long s);
	RNG(unsigned long init_key[], int key_length);
	~RNG() = default;

	// initializes mt[N] with a seed
	void seed(unsigned long s);

	void seed(unsigned long init_key[], int key_length);

	unsigned long randInt32();
	
	long randInt31();

	double randReal1();

	double randReal2();

	double randReal3();

	double randReal4();

	double randRes53();

	void demo();

private:
	const int N = 624;
	const int M = 397;
	const unsigned long MATRIX_A = 0x9908b0dfUL;
	const unsigned long LOWER_MASK = 0x7fffffffUL;

	int mti = 625;
	unsigned long mt[624];
};


}	// end namespace random
}	// end namespace minilab


class Demo
{
public:
	Demo();
};