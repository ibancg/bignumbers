#include <string.h>

#include "bignum.h"

// Addition. Simple implementation.
void add(BigNumber &A, BigNumber &B, BigNumber &C, bool sign) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive == B.isPositive) {

		// same sign case

		for (i = 0; i < N_DIGITS; i++) {
			r = carry + A.digits[i] + B.digits[i];
			carry = (r > 9) ? 1 : 0;
			C.digits[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;

	} else {

		// different sign case

		BigNumber* M; // higher module BN
		BigNumber* m; // lower module BN

		M = NULL;

		for (i = N_DIGITS - 1; i >= 0; i--) {

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
			memset(C.digits, 0, N_DIGITS * sizeof(bcd_t));
			C.isPositive = sign;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < N_DIGITS; i++) {
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
void sub(BigNumber &A, BigNumber &B, BigNumber &C, bool piz) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive != B.isPositive) {

		// different sign case

		for (i = 0; i < N_DIGITS; i++) {

			r = carry + A.digits[i] + B.digits[i];
			carry = (r > 9) ? 1 : 0;
			C.digits[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;
	} else {

		// same sign case

		BigNumber* M; // higher module BN
		BigNumber* m; // lower module BN

		M = NULL;

		for (i = N_DIGITS - 1; i >= 0; i--) {

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
			memset(C.digits, 0, N_DIGITS * sizeof(bcd_t));
			C.isPositive = piz;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < N_DIGITS; i++) {

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
