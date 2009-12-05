#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "bignum.h"

using namespace std;

// Computes the inverse.
// The problem is equivalent to find a root of the function f(x) = A*x - 1,
// so we can solve it by Newton's method doing the iteration x = x*(2 - A*x),
// that doesn't need any division.
void inv(BigNumber &A, BigNumber &B) {
	static BigNumber x1, x2;
	static BigNumber two("2");

# 	ifdef DEBUG
	cout << "INV";
	cout.flush();
#	endif

	B.isPositive = A.isPositive;
	memset(B.digits, 0, N_DIGITS * sizeof(bcd_t)); // cleaning the result

	int ipc = findFirstNonZeroDigitIndex(A);

	// if A has order n, 1/A has order -n, so we choose 10^(-n) as a starting
	// point.
	B.digits[N_FRAC_DIGITS - (ipc - N_FRAC_DIGITS) - 1] = 1;

	for (int k = 0;; k++) {

		mul(A, B, x1);
		sub(two, x1, x2);
		mul(B, x2, x1);

# 		ifdef DEBUG
		cout << '.';
		cout.flush();
#		endif

		if (equals(B, x1)) // loop until convergence
			break;

		copy(x1, B);
	}

# 	ifdef DEBUG
	cout << endl;
#	endif
}

#ifdef INVERSE_DIV_ALGORITHM

void div(BigNumber &A, BigNumber &B, BigNumber &C) {
	inv(B, C);
	mul(A, C, C);
}

#else

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
#endif
