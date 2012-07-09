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

private:
	void parse(const char *);

public:

	static long N_DIGITS; // number of digits in the format.
	static long N_FRAC_DIGITS; // number of fractional digits.

	bcd_t* digits;
	bool isPositive; // positive/!negative flag

	// Constructors.

	// Creates an empty bignumber
	BigNumber();
	BigNumber(const BigNumber&);

	// Creates a bignumber from a string, example: N = BigNumber("-1786.059e36");
	BigNumber(const char *);
	BigNumber(double);

	~BigNumber();

	// Assignation operator
	BigNumber& operator=(const BigNumber&);

	// sets to 0
	void clear();

	double toDouble() const;

	// Visualization. The parameter threshold configures the limit in number of
	// digits below which all the digits are shown. If the number of digits is
	// greater than the threshold, a format like 1274...0246.5162...2134 is
	// used. When this short notation is used, shortNotationDigits digits are
	// depicted in each group
	void show(std::ostream& ostream = std::cout, int threshold = 15000,
			int shortNotationDigits = 9) const;

	// Tests if two BNs are equal
	friend bool equals(const BigNumber& A, const BigNumber& B);

	// Compares two BNs and returns the number of coincident digits
	friend int compare(const BigNumber& A, const BigNumber& B);

	// First non-zero digit index.
	int firstNonZeroDigitIndex() const;

	// Operations.

	// Computes C = A + B. If the result is zero, the sign can be explicitly set.
	// A, B and C are overlappables.
	friend void add(const BigNumber& A, const BigNumber& B, BigNumber& C,
			bool sign = true);

	// Computes C = A - B. If the result is zero, the sign can be explicitly set.
	// A, B and C are overlappables.
	friend void sub(const BigNumber& A, const BigNumber& B, BigNumber& C,
			bool sign = true);

	// Computes C = A*B
	// A, B and C are overlappables
	friend void mul(const BigNumber& A, const BigNumber& B, BigNumber& C);

	// Computes the inverse of A, B = 1/A
	friend void inv(const BigNumber& A, BigNumber& B);

	// Computes C = A/B
	friend void div(const BigNumber& A, const BigNumber& B, BigNumber& C);

	// Computes the square root
	friend void sqrt(const BigNumber& A, BigNumber& B);

	// Computes the quartic root
	friend void sqrt4(const BigNumber& A, BigNumber& B);
};

#endif
