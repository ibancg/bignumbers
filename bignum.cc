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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// Generic functions implementation

// Constructs an empty BN.
BigNumber::BigNumber() {
}

// Constructs a BN by parsing the string representation
BigNumber::BigNumber(const char *s) {
	memset(digits, 0, N_DIGITS * sizeof(bcd_t));

	// if the first char is the minus sign, the number is negative
	isPositive = (*s != '-');

	if (!isPositive)
		s++;

	char *e = strchr(s, 'e');
	unsigned int ls = ((e) ? (e - s) : strlen(s));

	char *bp;
	int exponent = ((e) ? strtol(&e[1], &bp, 0) : 0);

	char *dot = strchr(s, '.');
	unsigned int dot_index = ((dot) ? (dot - s) : ls); // dot index.
	unsigned int i;

	for (i = 0; i < dot_index; i++)
		digits[N_FRAC_DIGITS + dot_index - i - 1 + exponent] = s[i] - 48;

	if (dot_index != ls)
		for (i = dot_index + 1; i < ls; i++)
			digits[N_FRAC_DIGITS - (i - dot_index) + exponent] = s[i] - 48;
}

void BigNumber::show() {
	int i, j;
	unsigned short int nc = 0;
	bool z = true;

	if (!isPositive)
		printf("-");

	for (i = N_DIGITS - 1; i >= N_FRAC_DIGITS; i--) {

		if (z && (digits[i]))
			z = false;
		if (!z) {
			printf("%c", digits[i] + 48);
			nc++;
		}
	}

	if (z) { // special case: zero.
		printf("0");
		nc++;
	}

	for (j = 0; j < N_FRAC_DIGITS; j++)
		if (digits[j])
			break;

	// decimal part.
	if (j != N_FRAC_DIGITS)
		printf(".");
	for (i = N_FRAC_DIGITS - 1; i >= j; i--) {
		printf("%c", digits[i] + 48);
		nc++;
	}

	printf("::(%u digits)\n", nc);
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
