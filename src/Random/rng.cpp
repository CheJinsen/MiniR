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

#include "Random/rng.h"


namespace minilab
{
namespace random
{


RNG::RNG()
{
	std::random_device rd;
	seed(rd());
}

RNG::RNG(unsigned long s)
{
	seed(s);
}

RNG::RNG(unsigned long init_key[], int key_length)
{
	seed(init_key, key_length);
}

void RNG::seed(unsigned long s)
{
	mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti-1] >> 30)) + mti); 
        mt[mti] &= 0xffffffffUL;
    }
}

void RNG::seed(unsigned long init_key[], int key_length)
{
    seed(19650218UL);
    int i = 1, j = 0;
    int k = N > key_length ? N : key_length;
    
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) + init_key[j] + j;
        mt[i] &= 0xffffffffUL;
        i++; j++;
        
        if (i >= N) {
            mt[0] = mt[N-1];
            i = 1;
        }

        if (j >= key_length) {
            j = 0;
        }
    }

    for (k = N - 1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i;
        mt[i] &= 0xffffffffUL;
        i++;

        if (i>=N) {
            mt[0] = mt[N - 1];
            i = 1;
        }
    }

    mt[0] = 0x80000000UL;
}

unsigned long RNG::randInt32()
{
	unsigned long y;
    unsigned long mag01[2] = {0x0UL, MATRIX_A};

    if (mti >= N) {
        int kk;

        if (mti == N+1) {
            seed(5489UL);
    	}

        for (kk = 0; kk < N - M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }

        for (; kk < N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }

        y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

long RNG::randInt31()
{
	return static_cast<long>(randInt32() >> 1);
}

double RNG::randReal1()
{
	// divided by 2^32 - 1
	return randInt32() * (1.0 / 4294967295.0);
}

double RNG::randReal2()
{
	// divided by 2^32
	return randInt32() * (1.0 / 4294967296.0);
}

double RNG::randReal3()
{
	// divided by 2^32
	return (static_cast<double>(randInt32()) + 0.5) * (1.0 / 4294967296.0);
}

double RNG::randRes53()
{
	unsigned long a = randInt32() >> 5;
	unsigned long b = randInt32() >> 6;
	return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}

}	// end namespace random
}	// end namespace minilab