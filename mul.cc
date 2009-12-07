/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2009, Ibán Cereijo Graña <ibancg at gmail dot com>

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

// Multiplication algorithms
//
//          AAAA.AAAA  The multiplication result has the double of digits than
//		  * BBBB.BBBB  the operands. In a fixed-point system like this, as
//  -----------------  digits as assigned decimal figures will be lost (by the
//  CCCCCCCC.CCCCCCCC  least significant side). In the most signiticant side,
//					   other digits are also dropped, but this is not considered
//					   as a loss of information, since if any of them were
//					   different to 0, an overflow error may happened.

#include "bignum.h"

// Classical algorithm.
void mul(BigNumber &A, BigNumber &B, BigNumber &digits) {
	int i, j;
	unsigned long int r, c;
	static bcd_t R[N_DIGITS + N_FRAC_DIGITS];

	memset(R, 0, (N_DIGITS + N_FRAC_DIGITS) * sizeof(bcd_t));

	for (i = 0; i < N_DIGITS; i++) {
		for (j = 0, c = 0; (j < N_DIGITS) && ((j + i) < (N_DIGITS
				+ N_FRAC_DIGITS)); j++) {

			r = R[j + i] + c + A.digits[j] * B.digits[i];
			c = (r / 10);
			R[j + i] = (r - 10 * c);
		}

		// decimal adjustment
		for (; (j + i) < (N_DIGITS + N_FRAC_DIGITS); j++) {
			r = R[j + i] + c;
			c = (r / 10);
			R[j + i] = (r - 10 * c);
		}
	}

	digits.isPositive = !(A.isPositive ^ B.isPositive);
	memcpy(digits.digits, &R[N_FRAC_DIGITS], N_DIGITS * sizeof(bcd_t));
}

