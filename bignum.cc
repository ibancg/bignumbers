/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2010, Ib�n Cereijo Gra�a <ibancg at gmail dot com>

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
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// Generic functions implementation

// Constructs an empty BN.
BigNumber::BigNumber() {
	digits = 0;
}

// Constructs a BN by parsing the string representation
BigNumber::BigNumber(const char *s) {
	digits = new bcd_t[N_DIGITS];
	memset(digits, 0, N_DIGITS * sizeof(bcd_t));

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
		digits[N_FRAC_DIGITS + dot_index - i - 1 + exponent] = s[i] - 48;

	int index;
	if (dot_index != ls)
		for (i = dot_index + 1; i < ls; i++) {
			index = N_FRAC_DIGITS - (i - dot_index) + exponent;
			if (index > 0) {
				digits[index] = s[i] - 48;
			}
		}
}

// Constructs an empty BN.
BigNumber::~BigNumber() {
	if (digits != 0) {
		delete[] digits;
	}
}

void BigNumber::show(int threshold, int shortNotationDigits) {
	int i, j;
	long int ni = 0; // number of integer digits
	long int nf = 0; // number of fractional digits
	bool z = true;

	if (!isPositive)
		printf("-");

	int firstNonZeroIndex = N_DIGITS;

	for (i = N_DIGITS - 1; i >= N_FRAC_DIGITS; i--) {

		if (z && (digits[i])) {
			firstNonZeroIndex = i;
			z = false;
		}
	}

	if (z) { // special case: zero.
		printf("0");
		ni = 0;
	} else {

		ni = firstNonZeroIndex - N_FRAC_DIGITS + 1;
		if ((ni > threshold) && (threshold > 2 * shortNotationDigits)) {
			for (i = firstNonZeroIndex; i > firstNonZeroIndex
					- shortNotationDigits; i--) {
				printf("%c", digits[i] + 48);
			}
			printf("...");
			for (i = N_FRAC_DIGITS + shortNotationDigits - 1; i
					>= N_FRAC_DIGITS; i--) {
				printf("%c", digits[i] + 48);
			}

		} else {
			for (i = firstNonZeroIndex; i >= N_FRAC_DIGITS; i--) {
				printf("%c", digits[i] + 48);
			}
		}

	}

	for (j = 0; j < N_FRAC_DIGITS; j++)
		if (digits[j])
			break;

	nf = N_FRAC_DIGITS - j;

	if (nf > 0) {
		// decimal part.
		printf(".");

		nf = N_FRAC_DIGITS - j;
		if ((nf > threshold) && (threshold > 2 * shortNotationDigits)) {
			for (i = N_FRAC_DIGITS - 1; i >= N_FRAC_DIGITS
					- shortNotationDigits; i--) {
				printf("%c", digits[i] + 48);
			}
			printf("...");
			for (i = j + shortNotationDigits - 1; i >= j; i--) {
				printf("%c", digits[i] + 48);
			}
		} else {
			for (i = N_FRAC_DIGITS - 1; i >= j; i--) {
				printf("%c", digits[i] + 48);
			}
		}
	}

	printf("::(%lu digits, %lu integer and %lu fractional)\n", ni + nf, ni, nf);
}

void copy(BigNumber &A, BigNumber &B) {
	memcpy(B.digits, A.digits, N_DIGITS * sizeof(bcd_t));
	B.isPositive = A.isPositive;
}

int findFirstNonZeroDigitIndex(BigNumber &A) {
	for (register int i = N_DIGITS - 1; i >= 0; i--)
		if (A.digits[i])
			return i;
	return -1; // special case: zero
}

bool equals(BigNumber &A, BigNumber &B) {
	for (register int i = 0; i < N_DIGITS; i++)
		if (A.digits[i] != B.digits[i])
			return false;
	return true;
}

int compare(BigNumber &A, BigNumber &B) {
	for (register int i = N_DIGITS - 1; i >= 0; i--)
		if (A.digits[i] != B.digits[i])
			return (N_DIGITS - i);

	return N_DIGITS;
}
