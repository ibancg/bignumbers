/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2010, Ibán Cereijo Graña <ibancg at gmail dot com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <sys/time.h>

#include "fft.h"
#include "bignum.h"

int main() {

	struct timeval t1, t2, t3;
	double elapsed_time;

	// initialize the fft library
	createPhaseFactors();

	// computes the Mersenne number 2^p - 1

	unsigned long int p = 3021377; // exponent (37th known Mersenne number)

	unsigned int nbits = 0; // number of bits of p
	unsigned long int p2 = p;
	while (p2 != 0) { // count the number of bits of p
		p2 >>= 1;
		nbits++;
	}

	unsigned int nmuls = 0; // number of needed muls

	BigNumber X = BigNumber("1");
	BigNumber AX = BigNumber("2");

	gettimeofday(&t1, NULL);
	printf("Evaluating the Mersenne number 2^%lu\n", p);
	printf("0%% completed\n");

	// the algorithm finds first the binary representation of p = (b0 b1 ... bn)
	// so p = b0*2^0 + b1*2^1 + ... + bn*2^n, that is,
	// 2^p = 2^(b0*2^0)*2^(b2*2^1)*...*2^(bn*2^n)
	// the 2^(2^i) quantity can be easily computed in a loop by multiplying an
	// acumulator by itself in each iteration (starting with 2)
	p2 = p;
	int i = 0;
	while (p2 != 0) { // starts process
		int bit = p2 & 1;
		p2 >>= 1;

		if (bit) {
			mulFFT(X, AX, X);
			nmuls++;
		}

		if (p == 0) {
			break;
		}

		mulFFT(AX, AX, AX);
		nmuls++;

		i++;
		printf("%d%% completed\n", (100 * i) / nbits);
	}

	BigNumber One = BigNumber("1");
	sub(X, One, X);

	gettimeofday(&t2, NULL);

	printf("result = ");
	X.show();
	timersub(&t2, &t1, &t3);
	elapsed_time = t3.tv_sec + 1e-6 * t3.tv_usec;
	printf("computation time: %lf seconds\n", elapsed_time);
	printf("%d multiplications and 1 substraction needed\n", nmuls);
}
