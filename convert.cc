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

// Floating point and BigNumber conversions.

// number of significant digits in floating point conversion.
#define N_SIGNIFICANT_DIGITS 20

void flt2BigNumber(flt_t x, BigNumber &X) {
	static char buffer[40];
	char *s = buffer;

	// prints the number in a string
	sprintf(s, "%20.20f", x);
	memset(X.digits, 0, N_DIGITS * sizeof(bcd_t));

	// and parses the string

	// if the string starts with a minus sign, the number is negative
	X.isPositive = (*s != '-');

	if (!X.isPositive)
		s++;

	char* e = strchr(s, 'e');
	unsigned int i;
	unsigned int ls = ((e) ? (e - s) : strlen(s));

	char* bp;
	int exponent = ((e) ? strtol(&e[1], &bp, 0) : 0);

	char* dot = strchr(s, '.');
	unsigned int dot_index = ((dot) ? (dot - s) : ls); // posición del punto.

	for (i = 0; i < dot_index; i++)
		X.digits[N_FRAC_DIGITS + dot_index - i - 1 + exponent] = s[i] - 48;

	if (dot_index != ls)
		for (i = dot_index + 1; i < ls; i++)
			X.digits[N_FRAC_DIGITS - (i - dot_index) + exponent] = s[i] - 48;
}

flt_t bigNumber2Flt(BigNumber &X) {
	register int i;
	double x;

	int n = findFirstNonZeroDigitIndex(X);

	if (n == -1) {
		return 0.0;
	}

	x = 0.0;
	i = n - N_SIGNIFICANT_DIGITS + 1;
	if (i < 0)
		i = 0;

	double exponent = pow(10, i - N_FRAC_DIGITS);

	for (; i <= n; i++) {
		x += X.digits[i] * exponent;
		exponent *= 10.0;
	}

	return x;
}
