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
#include <stdlib.h>
#include <algorithm>

#include "bignum.h"

using namespace std;

// Computes the inverse.
// The problem is equivalent to find a root of the function f(x) = 1/x - 1/A,
// so we can solve it by Newton's method doing the iteration x = x*(2 - A*x),
// that doesn't need any division.
void inv(const BigNumber &A, BigNumber &B) {
	static BigNumber x1, x2;
	static BigNumber two = 2.0;

# 	ifdef DEBUG
	cout << "INV";
	cout.flush();
#	endif

	double a;
	long int aexp;
	A.toDouble(a, aexp);
	B.fromDouble(1 / a, -aexp);

	//	int ipc = A.firstNonZeroDigitIndex();
	//
	//	B.clear(); // cleaning the result
	//	B.isPositive = A.isPositive;
	//
	//
	//	// TODO: improve first guess
	//	// if A has order n, 1/A has order -n, so we choose 10^(-n) as a starting
	//	// point.
	//	B.digits[BigNumber::N_FRAC_DIGITS - (ipc - BigNumber::N_FRAC_DIGITS) - 1] =
	//			1;

	for (int k = 0;; k++) {

		mul(A, B, x1);
		sub(two, x1, x2);
		mul(B, x2, x1);

# 		ifdef DEBUG
		cout << '.';
		cout.flush();
#		endif

		if (B == x1) // loop until convergence
			break;

		B = x1;
	}

# 	ifdef DEBUG
	cout << endl;
#	endif
}

// Computes the division by the classical Long Division Algorithm.
// A is overlappable with B or C.
void divLDA(const BigNumber &A, const BigNumber &B, BigNumber &C) {
	int i, j, n;
	int na = 0, nb = 0, nab;
	static std::vector<bcd_t> AA(
			BigNumber::N_DIGITS + BigNumber::N_FRAC_DIGITS);

	// AA is a shifted copy of A.
	copy(A.digits.begin(), A.digits.end(),
			AA.begin() + BigNumber::N_FRAC_DIGITS);
	fill(AA.begin(), AA.begin() + BigNumber::N_FRAC_DIGITS - 1, 0);

	for (i = BigNumber::N_DIGITS - 1; i >= 0; i--) {
		if ((!na) && (A.digits[i]))
			na = i;
		if ((!nb) && (B.digits[i]))
			nb = i;
	}

	// AA is shifted BigNumber::N_FRAC_DIGITS digits.
	na += BigNumber::N_FRAC_DIGITS;

	if (!nb && !B.digits[0]) {
		printf("ERROR: division by 0\n");
		exit(255);
	}

	C.clear();
	C.positive = !(A.positive ^ B.positive);

	if (nb > na)
		return;

	nab = na - nb;

	bool menor;
	char r, c;

	i = 1;
	for (i = 0; i <= nab; i++) {

		// Tests if AA is lower than a shifted version of B.
		for (n = 0;; n++) {

			menor = false;
			for (j = nb + 1; j >= 0; j--) {
				if (AA[j + nab - i] == B.digits[j])
					continue;
				if (AA[j + nab - i] < B.digits[j])
					menor = true;
				break;
			}
			if (menor)
				break;

			// Substracts a shifted version of B from AA.
			for (c = 0, j = 0; j <= nb + 1; j++) {
				r = AA[j + nab - i] - (B.digits[j] + c);
				c = (r < 0) ? 1 : 0;
				AA[j + nab - i] = (r + 10 * c); // r % 10
			}
		}

		C.digits[nab - i] = n;
	}
}

void div(const BigNumber &A, const BigNumber &B, BigNumber &C) {
#ifdef INVERSE_DIV_ALGORITHM
	inv(B, C);
	mul(A, C, C);
#else
	divLDA(A, B, C);
#endif
}

