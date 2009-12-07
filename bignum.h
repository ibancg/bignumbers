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

#ifndef __BIGNUMBER_H_
#define __BIGNUMBER_H_

// Arbitrary-precission fixed-point arithmetic.
//
// A BCD system with sign and modulus representation and decimal adjustment
// has been chosen, so each figure represents a decimal digit.

#include "config.h"

class BigNumber {

public:

	bcd_t digits[N_DIGITS];
	bool isPositive; // positive/!negative flag

	// Constructors.

	// Creates an empty bignumber
	BigNumber();

	// Creates a bignumber from a string, example: N = BigNumber("-1786.059e36");
	BigNumber(const char *);

	// Visualization.
	void show();

	// Conversions between bignumbers and floats
	friend void flt2BigNumber(flt_t flt, BigNumber& bn);
	friend flt_t bigNumber2Flt(BigNumber& bn);

	// Tests if two BNs are equal
	friend bool equals(BigNumber& A, BigNumber& B);

	// Compares two BNs and returns the number of coincident digits
	friend int compare(BigNumber& A, BigNumber& B);

	// First non-zero digit index.
	friend int findFirstNonZeroDigitIndex(BigNumber& A);

	// Copies a BN into another
	friend void copy(BigNumber& A, BigNumber& B);

	// Operations.

	// Computes C = A + B. If the result is zero, the sign can be explicitly set.
	// A, B and C are overlappables.
	friend void add(BigNumber& A, BigNumber& B, BigNumber& C, bool sign = true);

	// Computes C = A - B. If the result is zero, the sign can be explicitly set.
	// A, B and C are overlappables.
	friend void sub(BigNumber& A, BigNumber& B, BigNumber& C, bool sign = true);

	// Computes C = A*B
	// A, B and C are overlappables
	friend void mul(BigNumber& A, BigNumber& B, BigNumber& C);

};

#endif
