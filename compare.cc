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

#include "fft.h"
#include "bignum.h"
#include <stdlib.h>

#include <vector>

using namespace std;

void f_mul_fft(BigNumber& a, BigNumber& b, BigNumber& c) {
	mulFFT(a, b, c);
}

void f_mul_lma(BigNumber& a, BigNumber& b, BigNumber& c) {
	mulLMA(a, b, c);
}

void f_div_inv(BigNumber& a, BigNumber& b, BigNumber& c) {
	divInv(a, b, c);
}

void f_div_lda(BigNumber& a, BigNumber& b, BigNumber& c) {
	divLDA(a, b, c);
}

void f_inv(BigNumber& a, BigNumber& b, BigNumber&) {
	inv(a, b);
}

void f_sqrt_inv(BigNumber& a, BigNumber& b, BigNumber&) {
	sqrtInv(a, b);
}

void f_sqrt_noinv(BigNumber& a, BigNumber& b, BigNumber&) {
	sqrtNoInv(a, b);
}

void f_sqrt4_inv(BigNumber& a, BigNumber& b, BigNumber&) {
	sqrt4Inv(a, b);
}

void f_sqrt4_noinv(BigNumber& a, BigNumber& b, BigNumber&) {
	sqrt4NoInv(a, b);
}

void compareMethods() {

	std::string fileName = "output.txt";
	struct timeval t1, t2, t3;
	double elapsed_time;

	std::cout << "------- PI number computation -------" << std::endl;
	std::cout << "Computing first iterants ..." << std::endl;

	typedef void (func_t)(BigNumber& a, BigNumber& b, BigNumber&);

	vector<func_t*> fun;
	fun.push_back(&f_mul_fft);
	fun.push_back(&f_mul_lma);
	fun.push_back(&f_div_inv);
	fun.push_back(&f_div_lda);
	fun.push_back(&f_inv);
	fun.push_back(&f_sqrt_inv);
	fun.push_back(&f_sqrt_noinv);
	fun.push_back(&f_sqrt4_inv);
	fun.push_back(&f_sqrt4_noinv);
	int n = fun.size();
	vector<double> times(n);

	const int minExp = 5;
	const int maxExp = 14;
	const int maxRep = 4;

	const int innerReps[] = { 0, 0, 10000, 10000, 1000, 1000, 1000, 100, 100,
			10, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	const int maxExpPerFunction[] = { 30, 20, 30, 20, 30, 30, 25, 30, 25 };

	for (int exp = minExp; exp <= maxExp; exp++) {

		std::cout << "exp = " << exp << std::endl;

		const int innerRep = (exp < 20) ? innerReps[exp] : 1;

		BigNumber::N_DIGITS = (1 << exp);
		BigNumber::N_FRAC_DIGITS = BigNumber::N_DIGITS * 0.95;

		BigNumber a, b, c;

//		fill(times.begin(), times.end(), 0.0);

		for (int rep = 0; rep < maxRep; rep++) {
			for (long int i = 0; i < a.getNFracDigits() + 1; i++) {
				a[i] = rand() % 10;
				b[i] = rand() % 10;
			}

//			a.show();
//			b.show();

			for (int fi = 0; fi < n; fi++) {

				if (maxExpPerFunction[fi] >= exp) {
					gettimeofday(&t1, NULL);

					for (int ir = 0; ir < innerRep; ir++) {
						fun[fi](a, b, c);
					}

					gettimeofday(&t2, NULL);
					timersub(&t2, &t1, &t3);
					elapsed_time = t3.tv_sec + 1e-6 * t3.tv_usec;

					cout << "computation time: " << elapsed_time << " seconds"
							<< endl;
					times[fi] = elapsed_time / innerRep;
				} else {
					times[fi] = -1.0;
				}
			}

			std::ofstream file;
			file.open(fileName.c_str(), ios_base::app);
			file << a.getNFracDigits() + 1;
			for (unsigned int i = 0; i < times.size(); i++) {
				file << "   " << times[i];
			}
			file << "   " << endl;
			file.close();
		}

//		cout << "average times ";
//		for (unsigned int i = 0; i < times.size(); i++) {
//			cout << "   " << times[i]/maxRep;
//		}
//		cout << "   " << endl;

	}

}

