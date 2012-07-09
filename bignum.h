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

#ifndef __BIGNUMBER_H_
#define __BIGNUMBER_H_

// Arbitrary-precission fixed-point arithmetic.
//
// A BCD system with sign and modulus representation and decimal adjustment
// has been chosen, so each figure represents a decimal digit.

#include <ostream>
#include <iostream>
#include "config.h"

class BigNumber {

public:

	static long N_DIGITS; // number of digits in the format.
	static long N_FRAC_DIGITS; // number of fractional digits.

	bcd_t* digits;
	bool isPositive; // positive/!negative flag

	// Constructors.

	// Creates an empty bignumber
	BigNumber();

	// Creates a bignumber from a string, example: N = BigNumber("-1786.059e36");
	BigNumber(const char *);

	~BigNumber();

	// Assignation operator
	BigNumber& operator=(const BigNumber&);


	// Visualization. The parameter threshold configures the limit in number of
	// digits below which all the digits are shown. If the number of digits is
	// greater than the threshold, a format like 1274...0246.5162...2134 is
	// used. When this short notation is used, shortNotationDigits digits are
	// depicted in each group
	void show(std::ostream& ostream = std::cout, int threshold = 15000,
			int shortNotationDigits = 9);

	// Conversions between bignumbers and floats
	friend void flt2BigNumber(flt_t flt, BigNumber& bn);
	friend flt_t bigNumber2Flt(BigNumber& bn);

	// Tests if two BNs are equal
	friend bool equals(BigNumber& A, BigNumber& B);

	// Compares two BNs and returns the number of coincident digits
	friend int compare(BigNumber& A, BigNumber& B);

	// First non-zero digit index.
	friend int findFirstNonZeroDigitIndex(BigNumber& A);

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

	// Computes the inverse of A, B = 1/A
	friend void inv(BigNumber& A, BigNumber& B);

	// Computes C = A/B
	friend void div(BigNumber& A, BigNumber& B, BigNumber& C);

	// Computes the square root
	friend void sqrt(BigNumber& A, BigNumber& B);

	// Computes the quartic root
	friend void sqrt4(BigNumber& A, BigNumber& B);
};

#endif
