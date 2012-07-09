#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "bignum.h"

using namespace std;

// Square root computation by Newton's method. (F(x) = x^2 - 1/A).  This method
// converges to the inverse of the square root and doesn't need to compute any
// division in the time-critical loop.
// A and x NO overlappables.
void sqrtInv(const BigNumber &A, BigNumber &x) {
	static BigNumber xo = 0.0;
	static BigNumber _1p2 = 0.5;
	static BigNumber _3 = 3.0;

# ifdef DEBUG
	cout << "SQRT";
	cout.flush();
# endif

	if (!A.isPositive) {
		printf("ERROR: complex root.\n");
		exit(-1);
	}

	xo.clear();

	// if A has order n, we start the iteration at 10^(-n/2)
	xo.digits[BigNumber::N_FRAC_DIGITS
			- (A.firstNonZeroDigitIndex() - BigNumber::N_FRAC_DIGITS + 1) / 2] =
			1;

	for (int k = 0;; k++) {

		x = xo;
		mul(x, x, xo);
		mul(xo, A, xo); // A*x^2
		sub(xo, _3, xo, false); // (A*x^2 - 3)
		mul(_1p2, xo, xo); // 0.5*(A*x^2 - 3)
		mul(x, xo, xo); // 0.5*x*(A*x^2 - 3)
		xo.isPositive = !xo.isPositive;

#   ifdef DEBUG
		cout << '.';
		cout.flush();
#   endif

		if (compare(x, xo) >= (BigNumber::N_DIGITS - 2))
			break; // convergence.
	}

# ifdef DEBUG
	cout << endl;
# endif

	xo.isPositive = true; // rule out the negative solution

	// the method has converged to the inverse of the solution 1/sqrt(A). If
	// we multiply by A, we get A/sqrt(A) = sqrt(A).
	mul(A, xo, x);
}

// Square root computation by Newton's method. (F(x) = x^2 - A).
// A and x NO overlappables.
void sqrtNoInv(const BigNumber &A, BigNumber &x) {
	static BigNumber x2, Fx, DFx, xo;

# ifdef DEBUG
	cout << "SQRT";
	cout.flush();
# endif

	if (!A.isPositive) {
		printf("ERROR: complex root.\n");
		exit(-1);
	}

	xo.clear();

	// if A has order n, we start the iteration at 10^(-n/2)
	xo.digits[BigNumber::N_FRAC_DIGITS
			+ (A.firstNonZeroDigitIndex() - BigNumber::N_FRAC_DIGITS) / 2] = 1;

	for (;;) {

		x = xo;
		mul(x, x, x2);
		sub(x2, A, Fx, false); // F(x) = x^2 - A

		add(x, x, DFx); // F'(x) = 2*x
		div(Fx, DFx, x2); // F(x)/F'(x)
		sub(x, x2, xo); // x - F(x)/F'(x)

#   ifdef DEBUG
		cout << '.';
		cout.flush();
#   endif

		if (x == xo)
			break; // convergence.
	}

# ifdef DEBUG
	cout << endl;
# endif

	x.isPositive = true; // rule out the negative solution

	// if F(x) <= 0 -> x^2 <= A, we can exit. In other case, we need to
	// substract 1
	if (!Fx.isPositive)
		return; // F(x) <= 0

	// lets substract 1 from the result.
	char c = 1;

	// BCD adjustment.
	for (int i = 0;; i++) {
		x.digits[i] -= c;
		c = (((char) x.digits[i]) < 0) ? 1 : 0;
		if (!c)
			break; // if the carry is zero, we can exit.
		x.digits[i] += 10 * c;
	}
}

//---------------------------- QUARTIC ROOT ----------------------------

// Quartic root extraction by Newton's method. (F(x) = x^4 - 1/A). This method
// converges to the inverse of the quartic root and doesn't need to compute any
// division in the time-critical loop.
// A and x NO overlappables.
void sqrt4Inv(const BigNumber &A, BigNumber &x) {
	static BigNumber xo;
	static BigNumber _1p4 = 0.25;
	static BigNumber _5 = 5.0;

	if (!A.isPositive) {
		printf("ERROR: complex root.\n");
		exit(255);
	}

# ifdef DEBUG
	cout << "SQRT4";
	cout.flush();
# endif

	xo.clear();

	// if A has order n, we start the iteration at 10^(-n/2)
	xo.digits[BigNumber::N_FRAC_DIGITS
			- (A.firstNonZeroDigitIndex() - BigNumber::N_FRAC_DIGITS + 1) / 4] =
			1;

	for (int k = 0;; k++) {

		x = xo;
		mul(x, x, xo);
		mul(xo, xo, xo); // x^4
		mul(xo, A, xo); // A*x^4
		sub(xo, _5, xo, false); // (A*x^2 - 5)
		mul(_1p4, xo, xo); // 0.25*(A*x^4 - 5)
		mul(x, xo, xo); // 0.25*x*(A*x^4 - 5)
		xo.isPositive = !xo.isPositive;

#   ifdef DEBUG
		cout << '.';
		cout.flush();
#   endif

		if (compare(x, xo) >= (BigNumber::N_DIGITS - 2))
			break; // convergence.
	}

# ifdef DEBUG
	cout << endl;
# endif

	xo.isPositive = true; // rule out the negative solution.
	inv(xo, x);
}

// Quartic root extraction by Newton's method. (F(x) = x^4 - A).
// A and x NO overlappables.
void sqrt4NoInv(const BigNumber &A, BigNumber &x) {
	static BigNumber x1, x2, Fx, DFx, xo;

	if (!A.isPositive) {
		printf("ERROR: complex root.\n");
		exit(255);
	}

	xo.clear();

	// if A has order n, we start the iteration at 10^(-n/2)
	xo.digits[BigNumber::N_FRAC_DIGITS
			+ (A.firstNonZeroDigitIndex() - BigNumber::N_FRAC_DIGITS) / 4] = 1;

	for (;;) {

		x = xo;
		mul(x, x, x1); // x^2
		mul(x, x1, x1); // x^3
		mul(x, x1, x2); // x^4
		sub(x2, A, Fx, false); // F(x) = x^4 - A

		add(x1, x1, x1); // 2*x^3
		add(x1, x1, DFx); // F'(x) = 4*x^3
		div(Fx, DFx, x2); // F(x)/F'(x)
		sub(x, x2, xo); // x - F(x)/F'(x)

		if (x == xo)
			break; // convergence.
	}

	// Por si el método de Newton me converge a la solución negativa
	x.isPositive = true;

	// if F(x) <= 0 -> x^2 <= A, we can exit. In other case, we need to
	// substract 1
	if (!Fx.isPositive)
		return; // F(x) <= 0

	// lets substract 1 from the result.
	char c = 1;

	// BCD adjustment.
	for (int i = 0;; i++) {
		x.digits[i] -= c;
		c = (((char) x.digits[i]) < 0) ? 1 : 0;
		if (!c)
			break; // if the carry is zero, we can exit.
		x.digits[i] += 10 * c;
	}
}

void sqrt(const BigNumber &A, BigNumber &x) {
#ifdef INVERSE_NEWTON_SQRT_ALGORITHM
	sqrtInv(A, x);
#else
	sqrtNoInv(A, x);
#endif
}

void sqrt4(const BigNumber &A, BigNumber &x) {
#ifdef INVERSE_NEWTON_SQRT4_ALGORITHM
	sqrt4Inv(A, x);
#else
	sqrt4NoInv(A, x);
#endif
}
