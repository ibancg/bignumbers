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

	const char* x = "96171307032";
	unsigned int iter = 10;
	unsigned long int power = pow(2, iter);
	struct timeval t1, t2, t3;
	double elapsed_time;

	// initialize the fft library
	createPhaseFactors();

	BigNumber X = BigNumber(x);

	printf("Computing %s^%ld using the FFT multiplication algorithm ...\n", x,
			power);

	gettimeofday(&t1, NULL);

	// exponentiation in a simple way
	for (int i = 0; i < 10; i++) {
		mulFFT(X, X, X);
	}

	gettimeofday(&t2, NULL);

	printf("result = ");
	X.show();
	timersub(&t2, &t1, &t3);
	elapsed_time = t3.tv_sec + 1e-6*t3.tv_usec;
	printf("computation time: %lf seconds\n", elapsed_time);

	X = BigNumber(x);

	printf("\nComputing %s^%ld using the long multiplication algorithm ...\n",
			x, power);

	gettimeofday(&t1, NULL);

	// exponentiation in a simple way
	for (int i = 0; i < 10; i++) {
		mulLMA(X, X, X);
	}

	gettimeofday(&t2, NULL);

	printf("result = ");
	X.show();
	timersub(&t2, &t1, &t3);
	elapsed_time = t3.tv_sec + 1e-6*t3.tv_usec;
	printf("computation time: %lf seconds\n", elapsed_time);
}
