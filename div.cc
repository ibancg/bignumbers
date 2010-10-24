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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// Computes the division by the classical method.
// A is overlappable with B or C.
void div(BigNumber &A, BigNumber &B, BigNumber &C) {
	int i, j, n;
	int na = 0, nb = 0, nab;
	static bcd_t AA[N_DIGITS + N_FRAC_DIGITS];

	// AA is a shifted copy of A.
	memcpy(&AA[N_FRAC_DIGITS], A.digits, N_DIGITS * sizeof(bcd_t));
	memset(AA, 0, N_FRAC_DIGITS * sizeof(bcd_t));

	for (i = N_DIGITS - 1; i >= 0; i--) {
		if ((!na) && (A.digits[i]))
			na = i;
		if ((!nb) && (B.digits[i]))
			nb = i;
	}

	// AA is shifted N_FRAC_DIGITS digits.
	na += N_FRAC_DIGITS;

	if (!nb && !B.digits[0]) {
		printf("ERROR: division by 0\n");
		exit(255);
	}

	memset(C.digits, 0, N_DIGITS * sizeof(bcd_t));
	C.isPositive = !(A.isPositive ^ B.isPositive);

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
