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

#include <fstream>
#include <sys/time.h>
#include <getopt.h>
#include <stdlib.h>

#include "fft.h"
#include "bignum.h"

int main(int argc, char *argv[]) {

	const char* const op_short = "o";
	const struct option op_long[] = { //
			{ "output", 1, NULL, 'o' }, //
					{ "help", 0, NULL, 'h' } //
			};

	int ich;
	std::string fileName = "output.txt";
	char* output_file_name = 0x0;
	while ((ich = getopt_long(argc, argv, op_short, op_long, NULL)) != EOF) {
		switch (ich) {
		case 'o':
			output_file_name = optarg;
			if (output_file_name) {
				std::cout << "output_file_name = " << output_file_name
						<< std::endl;
				fileName = output_file_name;
			}
			break;
		case 'h':
			printf("Usage: bignum [-o output_file]");
			printf("\nMersenne prime number evaluation.\n");
			printf("Example: bignum -o output.txt\n\n");
			printf("Options:\n");
			printf("  -o, --output=FILE  write the solution to FILE\n");
			printf("  -h, --help         display this help and exit\n");
			printf("\n");
			exit(0);
		default:
			printf("Usage: bignum [-o output_file]\n");
			printf("Try bignum --help for more information.\n");
			exit(0);
		}
	}

	struct timeval t1, t2, t3;
	double elapsed_time;

	// initialize the fft library
	createPhaseFactors();

	// computes the Mersenne number 2^p - 1

	// exponent
	//unsigned long int p = 43112609; // largest known Mersenne number
	//unsigned long int p = 3021377; // 37th known Mersenne number
	unsigned long int p = 246; // 37th known Mersenne number

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
	std::cout << "Evaluating the Mersenne number 2^" << p << " - 1"
			<< std::endl;
	std::cout << "0% completed" << std::endl;

	// the algorithm finds first the binary representation of p = (bn ... b1 b0)
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

		mulFFT(AX, AX, AX); // 2^(2^i)
		nmuls++;

		i++;
		std::cout << ((100 * i) / nbits) << "% completed" << std::endl;
	}

	BigNumber One = BigNumber("1");
	sub(X, One, X);

	gettimeofday(&t2, NULL);

	std::cout << "result = ";
	X.show();
	timersub(&t2, &t1, &t3);
	elapsed_time = t3.tv_sec + 1e-6 * t3.tv_usec;
	std::cout << "computation time: " << elapsed_time << " seconds"
			<< std::endl;
	std::cout << nmuls << " multiplications and 1 subtraction needed"
			<< std::endl;
	std::cout << "dumping result to file '" << fileName << "' ..." << std::endl;
	std::ofstream file(fileName.c_str());
	X.show(file, 0);
	std::cout << "done." << std::endl;

	destroyPhaseFactors();
}
