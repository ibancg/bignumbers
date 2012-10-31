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
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "bignum.h"

long BigNumber::N_DIGITS = 32; // number of digits in the format.
long BigNumber::N_FRAC_DIGITS = 16; // number of fractional digits.

// Constructs an empty BN.
BigNumber::BigNumber(long int nDigits, long int nFracDigits) :
		digits(nDigits)  {
	fill(digits.begin(), digits.end(), 0);
	positive = true;
	this->nDigits = nDigits;
	this->nFracDigits = nFracDigits;
}

BigNumber::BigNumber(const BigNumber& b) {
	*this = b;
}

// Constructs a BN by parsing the string representation
void BigNumber::parse(const char *s) {

	clear();

	// if the first char is the minus sign, the number is negative
	positive = (*s != '-');

	if (!positive)
		s++;

	const char *e = strchr(s, 'e');
	unsigned int ls = ((e) ? (e - s) : strlen(s));

	char *bp;
	int exponent = ((e) ? strtol(&e[1], &bp, 0) : 0);

	const char *dot = strchr(s, '.');
	unsigned int dot_index = ((dot) ? (dot - s) : ls); // dot index.
	unsigned int i;

	for (i = 0; i < dot_index; i++)
		digits[nFracDigits + dot_index - i - 1 + exponent] = s[i]
				- 48;

	int index;
	if (dot_index != ls)
		for (i = dot_index + 1; i < ls; i++) {
			index = nFracDigits - (i - dot_index) + exponent;
			if (index > 0) {
				digits[index] = s[i] - 48;
			}
		}
}

BigNumber::BigNumber(const char *s, long int nDigits, long int nFracDigits) :
		digits(nDigits) {
	this->nDigits = nDigits;
	this->nFracDigits = nFracDigits;
	parse(s);
}

BigNumber::BigNumber(double b, long int nDigits, long int nFracDigits) :
		digits(nDigits) {
	this->nDigits = nDigits;
	this->nFracDigits = nFracDigits;
	fromDouble(b);
}

// Constructs an empty BN.
BigNumber::~BigNumber() {
}

BigNumber& BigNumber::operator=(const BigNumber& a) {
	digits = a.digits;
	positive = a.positive;
	nDigits = a.nDigits;
	nFracDigits = a.nFracDigits;

	return *this;
}

void BigNumber::clear() {
	fill(digits.begin(), digits.end(), 0);
	positive = true;
}

void BigNumber::resize(long int nDigits, long int nFracDigits) {
	digits.resize(nDigits);
	this->nDigits = nDigits;
	this->nFracDigits = nFracDigits;
}

void BigNumber::resize(const BigNumber& A) {
	resize(A.nDigits, A.nFracDigits);
}


double BigNumber::toDouble() const {
	double x;
	long int xexp;
	toDouble(x, xexp);
	return x * pow10(xexp);
}

void BigNumber::toDouble(double& x, long int& exp) const {
	register int i;

	int n = firstNonZeroDigitIndex();

	if (n == -1) {
		x = 0.0;
		exp = 0;
	} else {
		x = 0.0;
		i = n;
		static const int N_SIGNIFICANT_DIGITS = 20;
		int j = std::max<int>(n - N_SIGNIFICANT_DIGITS, 0);
		double exponent = 1.0;
		for (; i >= j; i--) {
			x += digits[i] * exponent;
			exponent *= 0.1;
		}
		if (!positive) {
			x = -x;
		}
		exp = n - nFracDigits;
	}
}

void BigNumber::fromDouble(double b) {
	// TODO: do not use strings
	static char buffer[40];
	// prints the number in a string
	sprintf(buffer, "%20.20f", b);
	parse(buffer);
}

void BigNumber::fromDouble(double b, long int exp) {
	// TODO: do not use strings
	static char buffer[30];
	// prints the number in a string
	sprintf(buffer, "%0.20fe%li", b, exp);
	parse(buffer);
}

bool BigNumber::isPositive_() const {
	return positive;
}

void BigNumber::show(std::ostream& ostream, int threshold,
		int shortNotationDigits) const {
	int i, j;
	long int ni = 0; // number of integer digits
	long int nf = 0; // number of fractional digits
	bool z = true;

	if (!positive) {
		ostream << '-';
	}

	int firstNonZeroIndex = nDigits;

	for (i = nDigits - 1; i >= nFracDigits; i--) {

		if (z && (digits[i])) {
			firstNonZeroIndex = i;
			z = false;
		}
	}

	if (z) { // special case: zero.
		ostream << '0';
		ni = 0;
	} else {

		ni = firstNonZeroIndex - nFracDigits + 1;
		if ((threshold > 0) && (ni > threshold)
				&& (threshold > 2 * shortNotationDigits)) {
			for (i = firstNonZeroIndex;
					i > firstNonZeroIndex - shortNotationDigits; i--) {
				ostream << (char) (digits[i] + 48);
			}
			ostream << "...";
			for (i = nFracDigits + shortNotationDigits - 1; i >= nFracDigits;
					i--) {
				ostream << (char) (digits[i] + 48);
			}

		} else {
			for (i = firstNonZeroIndex; i >= nFracDigits; i--) {
				ostream << (char) (digits[i] + 48);
			}
		}

	}

	for (j = 0; j < nFracDigits; j++)
		if (digits[j])
			break;

	nf = nFracDigits - j;

	if (nf > 0) {
		// decimal part.
		ostream << '.';

		nf = nFracDigits - j;
		if ((nf > threshold) && (threshold > 2 * shortNotationDigits)) {
			for (i = nFracDigits - 1; i >= nFracDigits - shortNotationDigits;
					i--) {
				ostream << (char) (digits[i] + 48);
			}
			ostream << "...";
			for (i = j + shortNotationDigits - 1; i >= j; i--) {
				ostream << (char) (digits[i] + 48);
			}
		} else {
			for (i = nFracDigits - 1; i >= j; i--) {
				ostream << (char) (digits[i] + 48);
			}
		}
	}

	ostream << "::(" << (ni + nf) << " digits, " << ni << " integer and " << nf
			<< " fractional)" << std::endl;
}

int BigNumber::firstNonZeroDigitIndex() const {
	for (register int i = nDigits - 1; i >= 0; i--)
		if (digits[i])
			return i;
	return -1; // special case: zero
}

bool operator==(const BigNumber &A, const BigNumber &B) {
	return A.digits == B.digits;
}

int matchingDigits(const BigNumber &A, const BigNumber &B) {
	long int i;
	for (i = A.nDigits - 1; i >= 0; i--)
		if (A.digits[i] != B.digits[i]) {
			return (A.nDigits - i);
		}
	return A.nDigits;
}

bool matchDimensions(const BigNumber& A, const BigNumber& B) {
	return (A.nDigits == B.nDigits) && (A.nFracDigits == B.nFracDigits);
}


bcd_t BigNumber::operator()(int i) const {
	return digits[nFracDigits + i];
}

bcd_t& BigNumber::operator()(int i) {
	return digits[nFracDigits + i];
}

bcd_t BigNumber::operator[](int i) const {
	return digits[i];
}

bcd_t& BigNumber::operator[](int i) {
	return digits[i];
}

long int BigNumber::getNDigits() const {
	return nDigits;
}

long int BigNumber::getNFracDigits() const {
	return nFracDigits;
}

long int BigNumber::getNIntDigits() const {
	return nDigits - nFracDigits;
}
