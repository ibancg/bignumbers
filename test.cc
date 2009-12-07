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

#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include "fft.h"
#include "bignum.h"

int main() {

	BigNumber X, Y, Z;

	createPhaseFactors();

	BigNumber X1, X2, X3, AX1, AX2, AX3, AX;

	printf("testing if 1841^12 + 1782^12 = 1922^12\t\n");

	X1 = BigNumber("1841");
	X2 = BigNumber("1782");
	X3 = BigNumber("1922");

	copy(X1, AX1);
	copy(X2, AX2);
	copy(X3, AX3);

	for (int i = 1; i < 12; i++) {

		mul(AX1, X1, AX);
		copy(AX, AX1);
		mul(AX2, X2, AX);
		copy(AX, AX2);
		mul(AX3, X3, AX);
		copy(AX, AX3);
	}

	printf("1841^12 = \t\t");
	AX1.show();
	printf("1782^12 =\t\t");
	AX2.show();
	add(AX1, AX2, AX);
	printf("1841^12 + 1782^12 =\t");
	AX.show();
	printf("1922^12 =\t\t");
	AX3.show();
	printf("%s\n", equals(AX, AX3) ? "true!" : "false!");
}
