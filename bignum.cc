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
BigNumber::BigNumber() :
		digits(BigNumber::N_DIGITS) {
	fill(digits.begin(), digits.end(), 0);
	isPositive = true;
}

BigNumber::BigNumber(const BigNumber& b) :
		digits(BigNumber::N_DIGITS) {
	*this = b;
}

// Constructs a BN by parsing the string representation
void BigNumber::parse(const char *s) {
	clear();

	// if the first char is the minus sign, the number is negative
	isPositive = (*s != '-');

	if (!isPositive)
		s++;

	const char *e = strchr(s, 'e');
	unsigned int ls = ((e) ? (e - s) : strlen(s));

	char *bp;
	int exponent = ((e) ? strtol(&e[1], &bp, 0) : 0);

	const char *dot = strchr(s, '.');
	unsigned int dot_index = ((dot) ? (dot - s) : ls); // dot index.
	unsigned int i;

	for (i = 0; i < dot_index; i++)
		digits[BigNumber::N_FRAC_DIGITS + dot_index - i - 1 + exponent] = s[i]
				- 48;

	int index;
	if (dot_index != ls)
		for (i = dot_index + 1; i < ls; i++) {
			index = BigNumber::N_FRAC_DIGITS - (i - dot_index) + exponent;
			if (index > 0) {
				digits[index] = s[i] - 48;
			}
		}
}

BigNumber::BigNumber(const char *s) :
		digits(BigNumber::N_DIGITS) {
	parse(s);
}

BigNumber::BigNumber(double b) :
		digits(BigNumber::N_DIGITS) {
	fromDouble(b);
}

// Constructs an empty BN.
BigNumber::~BigNumber() {
}

BigNumber& BigNumber::operator=(const BigNumber& a) {
	digits = a.digits;
	isPositive = a.isPositive;

	return *this;
}

void BigNumber::clear() {
	fill(digits.begin(), digits.end(), 0);
	isPositive = true;
}

double BigNumber::toDouble() const {
	register int i;
	double x;

	int n = firstNonZeroDigitIndex();

	if (n == -1) {
		return 0.0;
	}

	x = 0.0;
	static const int N_SIGNIFICANT_DIGITS = 20;
	i = n - N_SIGNIFICANT_DIGITS + 1;
	if (i < 0)
		i = 0;

	double exponent = pow(10, i - BigNumber::N_FRAC_DIGITS);

	for (; i <= n; i++) {
		x += digits[i] * exponent;
		exponent *= 10.0;
	}

	return x;
}

void BigNumber::fromDouble(double b) {
	static char buffer[40];
	// prints the number in a string
	sprintf(buffer, "%20.20f", b);
	parse(buffer);
}


void BigNumber::show(std::ostream& ostream, int threshold,
		int shortNotationDigits) const {
	int i, j;
	long int ni = 0; // number of integer digits
	long int nf = 0; // number of fractional digits
	bool z = true;

	if (!isPositive) {
		ostream << '-';
	}

	int firstNonZeroIndex = BigNumber::N_DIGITS;

	for (i = BigNumber::N_DIGITS - 1; i >= BigNumber::N_FRAC_DIGITS; i--) {

		if (z && (digits[i])) {
			firstNonZeroIndex = i;
			z = false;
		}
	}

	if (z) { // special case: zero.
		ostream << '0';
		ni = 0;
	} else {

		ni = firstNonZeroIndex - BigNumber::N_FRAC_DIGITS + 1;
		if ((threshold > 0) && (ni > threshold)
				&& (threshold > 2 * shortNotationDigits)) {
			for (i = firstNonZeroIndex;
					i > firstNonZeroIndex - shortNotationDigits; i--) {
				ostream << (char) (digits[i] + 48);
			}
			ostream << "...";
			for (i = BigNumber::N_FRAC_DIGITS + shortNotationDigits - 1;
					i >= BigNumber::N_FRAC_DIGITS; i--) {
				ostream << (char) (digits[i] + 48);
			}

		} else {
			for (i = firstNonZeroIndex; i >= BigNumber::N_FRAC_DIGITS; i--) {
				ostream << (char) (digits[i] + 48);
			}
		}

	}

	for (j = 0; j < BigNumber::N_FRAC_DIGITS; j++)
		if (digits[j])
			break;

	nf = BigNumber::N_FRAC_DIGITS - j;

	if (nf > 0) {
		// decimal part.
		ostream << '.';

		nf = BigNumber::N_FRAC_DIGITS - j;
		if ((nf > threshold) && (threshold > 2 * shortNotationDigits)) {
			for (i = BigNumber::N_FRAC_DIGITS - 1;
					i >= BigNumber::N_FRAC_DIGITS - shortNotationDigits; i--) {
				ostream << (char) (digits[i] + 48);
			}
			ostream << "...";
			for (i = j + shortNotationDigits - 1; i >= j; i--) {
				ostream << (char) (digits[i] + 48);
			}
		} else {
			for (i = BigNumber::N_FRAC_DIGITS - 1; i >= j; i--) {
				ostream << (char) (digits[i] + 48);
			}
		}
	}

	ostream << "::(" << (ni + nf) << " digits, " << ni << " integer and " << nf
			<< " fractional)" << std::endl;
}

int BigNumber::firstNonZeroDigitIndex() const {
	for (register int i = BigNumber::N_DIGITS - 1; i >= 0; i--)
		if (digits[i])
			return i;
	return -1; // special case: zero
}

bool operator==(const BigNumber &A, const BigNumber &B) {
	return A.digits == B.digits;
}

int compare(const BigNumber &A, const BigNumber &B) {
	for (register int i = BigNumber::N_DIGITS - 1; i >= 0; i--)
		if (A.digits[i] != B.digits[i])
			return (BigNumber::N_DIGITS - i);

	return BigNumber::N_DIGITS;
}
