/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2010, Ib�n Cereijo Gra�a <ibancg at gmail dot com>

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

#include "fft.h"
#include "bignum.h"

/**
 * Estimation of PI using a 4th order Borwein's algorithm (1985)
 */

void computePi() {

	std::string fileName = "output.txt";
	struct timeval t1, t2, t3;
	double elapsed_time;

	gettimeofday(&t1, NULL);

	std::cout << "------- PI number computation -------" << std::endl;
	std::cout << "Computing first iterants ..." << std::endl;

	BigNumber::N_DIGITS = (1 << 20);
	BigNumber::N_FRAC_DIGITS = BigNumber::N_DIGITS * 0.99;

	BigNumber UNO = 1.0;
	BigNumber ACP2 = 2.0; // powers of 2 accumulator.
	BigNumber y, x1, x2, x3, pio;
	bool stop;
	int i, j;

	BigNumber pi = 0.0;

	sqrt(ACP2, x2);
	sub(x2, UNO, y); // y0 = sqrt(2) - 1

	add(x2, x2, x2); // 2*sqrt(2)
	add(x2, x2, x2); // 4*sqrt(2)

	BigNumber a = 6.0;
	sub(a, x2, a);  // a0 = 6 - 4*sqrt(2)
	std::cout << "...OK" << std::endl;

	for (i = 0;; i++) {

		inv(a, pio);

		j = matchingDigits(pio, pi);
		stop = (j == pi.getNDigits());

		std::cout << "iteration " << (i + 1) << " : "
				<< j - (pi.getNIntDigits())
				<< " decimals found" << std::endl;

		if (stop)
			break;
		pi = pio;

		mul(y, y, x1);   // y^2
		mul(x1, x1, x2); // y^4
		sub(UNO, x2, x1); // 1 - y^4

#define USESQRT4

#ifdef	USESQRT4
		sqrt4(x1, x2); // (1 - y^4)^(1/4);
		sub(UNO, x2, x1); // (1 - (1 - y^4)^(1/4))
		add(UNO, x2, x3);  // (1 + (1 - y^4)^(1/4))
		div(x1, x3, y);     // (1 - (1 - y^4)^(1/4))/(1 + (1 - y^4)^(1/4))
#else
				sqrt(x1, x2); // (1 - y^4)^(1/2);
				sqrt(x2, x1);// (1 - y^4)^(1/4);
				sub(UNO, x1, x2);// (1 - (1 - y^4)^(1/4))
				add(UNO, x1, x3);// (1 + (1 - y^4)^(1/4))
				div(x2, x3, y);// (1 - (1 - y^4)^(1/4))/(1 + (1 - y^4)^(1/4))
#endif

		add(y, UNO, x1); // (1 + y)
		mul(x1, x1, x2);  // (1 + y)^2
		mul(x2, x2, x3);  // (1 + y)^4
		mul(x3, a, x2);   // (1 + y)^4*a

		add(ACP2, ACP2, ACP2); // 2*ACP2
		add(ACP2, ACP2, ACP2); // 4*ACP2

		mul(y, y, x3); // y^2
		add(x1, x3, x1); // (1 + y + y^2)
		mul(y, x1, x3);     // y*(1 + y + y^2)
		mul(ACP2, x3, x1); // 2^(2*i + 1)*y*(1 + y + y^2)
		sub(x2, x1, a); // a*(1 + y)^4 - 2^(2*i + 1)*y*(1 + y + y^2)

	}

	gettimeofday(&t2, NULL);

	std::cout << "result = ";
	pi.show();
	std::cout << i << " iterations needed to find " << pi.getNFracDigits()
			<< " decimal digits of PI." << std::endl;
	timersub(&t2, &t1, &t3);
	elapsed_time = t3.tv_sec + 1e-6 * t3.tv_usec;
	std::cout << "computation time: " << elapsed_time << " seconds"
			<< std::endl;
	std::cout << "dumping result to file '" << fileName << "' ..." << std::endl;
	std::ofstream file(fileName.c_str());
	file << "PI ~= ";
	pi.show(file, 0);
	std::cout << "done." << std::endl;

}
