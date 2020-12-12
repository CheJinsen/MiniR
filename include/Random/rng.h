/* 
   Random generator for MT19937.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Copyright (C) 2005, Mutsuo Saito
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

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
	const unsigned long UPPER_MASK = 0x80000000UL;
	const unsigned long LOWER_MASK = 0x7fffffffUL;

	int mti = 625;
	unsigned long mt[624];
};


}	// end namespace random
}	// end namespace minilab