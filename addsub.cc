#include <string.h>

#include "bignum.h"

// Addition. Simple implementation.
void add(BigNumber &A, BigNumber &B, BigNumber &C, bool sign) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive == B.isPositive) {

		// same sign case

		for (i = 0; i < NCIF; i++) {
			r = carry + A.C[i] + B.C[i];
			carry = (r > 9) ? 1 : 0;
			C.C[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;

	} else {

		// different sign case

		BigNumber* M; // higher module BN
		BigNumber* m; // lower module BN

		M = NULL;

		for (i = NCIF - 1; i >= 0; i--) {

			if (A.C[i] == B.C[i])
				continue;

			if (A.C[i] > B.C[i]) {
				M = &A;
				m = &B;
			} else {
				M = &B;
				m = &A;
			}
			break;
		}

		if (!M) { // the both numbers have the same module, so the result is 0
			memset(C.C, 0, NCIF * sizeof(TBC));
			C.isPositive = sign;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < NCIF; i++) {
			r = M->C[i] - (m->C[i] + carry);
			carry = (r < 0) ? 1 : 0;
			C.C[i] = (r + 10 * carry);
		}

		// if the number with higher module is positive, then the result is also
		// positive.
		C.isPositive = ((A.isPositive) && (M == &A)) || ((B.isPositive) && (M
				== &B));
	}
}

// Substraction. SImple implementation.
void sub(BigNumber &A, BigNumber &B, BigNumber &C, bool piz) {
	register int i;
	char r;
	char carry = 0;

	if (A.isPositive != B.isPositive) {

		// different sign case

		for (i = 0; i < NCIF; i++) {

			r = carry + A.C[i] + B.C[i];
			carry = (r > 9) ? 1 : 0;
			C.C[i] = (r - 10 * carry); // r % 10
		}

		C.isPositive = A.isPositive;
	} else {

		// same sign case

		BigNumber* M; // higher module BN
		BigNumber* m; // lower module BN

		M = NULL;

		for (i = NCIF - 1; i >= 0; i--) {

			if (A.C[i] == B.C[i])
				continue;

			if (A.C[i] > B.C[i]) {
				M = &A;
				m = &B;
			} else {
				M = &B;
				m = &A;
			}
			break;
		}

		if (!M) { // the both numbers have the same module, so the result is 0
			memset(C.C, 0, NCIF * sizeof(TBC));
			C.isPositive = piz;
			return;
		}

		// substracts the lower module number from the higher module one
		for (i = 0; i < NCIF; i++) {

			r = M->C[i] - (m->C[i] + carry);
			carry = (r < 0) ? 1 : 0;
			C.C[i] = (r + 10 * carry);
		}

		// if the number with higher module is positive, then the result is also
		// positive
		C.isPositive = ((A.isPositive) && (M == &A)) || ((!B.isPositive) && (M
				== &B));
	}
}
