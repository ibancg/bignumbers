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

#include "bignum.h"

// Addition. Simple implementation.
void add(const BigNumber &A, const BigNumber &B, BigNumber &C, bool sign) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive == B.isPositive) {

		// same sign case

		for (i = 0; i < BigNumber::N_DIGITS; i++) {
			r = carry + A.digits[i] + B.digits[i];
			carry = (r > 9) ? 1 : 0;
			C.digits[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;

	} else {

		// different sign case

		const BigNumber* M = 0x0; // higher module BN
		const BigNumber* m = 0x0; // lower module BN

		for (i = BigNumber::N_DIGITS - 1; i >= 0; i--) {

			if (A.digits[i] == B.digits[i])
				continue;

			if (A.digits[i] > B.digits[i]) {
				M = &A;
				m = &B;
			} else {
				M = &B;
				m = &A;
			}
			break;
		}

		if (!M) { //  both numbers have the same module, so the result is 0
			C.clear();
			C.isPositive = sign;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < BigNumber::N_DIGITS; i++) {
			r = M->digits[i] - (m->digits[i] + carry);
			carry = (r < 0) ? 1 : 0;
			C.digits[i] = (r + 10 * carry);
		}

		// if the number with higher module is positive, then the result is also
		// positive.
		C.isPositive = ((A.isPositive) && (M == &A)) || ((B.isPositive) && (M
				== &B));
	}
}

// Substraction. Simple implementation.
void sub(const BigNumber &A, const BigNumber &B, BigNumber &C, bool piz) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive != B.isPositive) {

		// different sign case

		for (i = 0; i < BigNumber::N_DIGITS; i++) {

			r = carry + A.digits[i] + B.digits[i];
			carry = (r > 9) ? 1 : 0;
			C.digits[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;
	} else {

		// same sign case

		const BigNumber* M = 0x0; // higher module BN
		const BigNumber* m = 0x0; // lower module BN

		for (i = BigNumber::N_DIGITS - 1; i >= 0; i--) {

			if (A.digits[i] == B.digits[i])
				continue;

			if (A.digits[i] > B.digits[i]) {
				M = &A;
				m = &B;
			} else {
				M = &B;
				m = &A;
			}
			break;
		}

		if (!M) { // both numbers have the same module, so the result is 0
			C.clear();
			C.isPositive = piz;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < BigNumber::N_DIGITS; i++) {

			r = M->digits[i] - (m->digits[i] + carry);
			carry = (r < 0) ? 1 : 0;
			C.digits[i] = (r + 10 * carry);
		}

		// if the number with higher module is positive, then the result is also
		// positive
		C.isPositive = ((A.isPositive) && (M == &A)) || ((!B.isPositive) && (M
				== &B));
	}
}
