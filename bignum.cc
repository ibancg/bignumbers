#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// Generic functions implementation

// Constructs an empty BN.
BigNumber::BigNumber() {
}


// Constructs a BN by parsing the string representation
BigNumber::BigNumber(char *s) {
	memset(C, 0, NCIF * sizeof(TBC));

	// if the first char is the minus sign, the number is negative
	isPositive = (*s != '-');

	if (!isPositive)
		s++;

	char *e = strchr(s, 'e');
	unsigned int ls = ((e) ? (e - s) : strlen(s));

	char *bp;
	int exponent = ((e) ? strtol(&e[1], &bp, 0) : 0);

	char *dot = strchr(s, '.');
	int pp = ((dot) ? (dot - s) : ls); // posición del punto.
	int i;

	for (i = 0; i < pp; i++)
		C[NFRC + pp - i - 1 + exponent] = s[i] - 48;

	if (pp != ls)
		for (i = pp + 1; i < ls; i++)
			C[NFRC - (i - pp) + exponent] = s[i] - 48;
}


void BigNumber::show() {
	int i, j;
	unsigned short int nc = 0;
	unsigned long int c, ch;
	bool z = true;

	if (!isPositive)
		printf("-");

	for (i = NCIF - 1; i >= NFRC; i--) {

		if (z && (C[i]))
			z = false;
		if (!z) {
			printf("%c", C[i] + 48);
			nc++;
		}
	}

	if (z) { // special case: zero.
		printf("0");
		nc++;
	}

	for (j = 0; j < NFRC; j++)
		if (C[j])
			break;

	// decimal part.
	if (j != NFRC)
		printf(".");
	for (i = NFRC - 1; i >= j; i--) {
		printf("%c", C[i] + 48);
		nc++;
	}

	printf("::(%u digits)\n", nc);
}


void shl(BigNumber &A, BigNumber &B, int d) {
	if (d >= 0) {
		memcpy(B.C + d, A.C, NCIF - d * sizeof(TBC));
		memset(B.C, 0, d * sizeof(TBC));
	} else {
		memcpy(B.C, A.C - d, NCIF + d * sizeof(TBC));
		memset(&B.C[NCIF + d], 0, -d * sizeof(TBC));
	}

	B.isPositive = A.isPositive;
}


void shr(BigNumber &A, BigNumber &B, int d) {
	shl(A, B, -d);
}


void copy(BigNumber &A, BigNumber &B) {
	memcpy(B.C, A.C, NCIF * sizeof(TBC));
	B.isPositive = A.isPositive;
}


int findFirstNonZeroDigitIndex(BigNumber &A) {
	for (register int i = NCIF - 1; i >= 0; i--)
		if (A.C[i])
			return i;
	return -1; // special case: zero
}


bool equals(BigNumber &A, BigNumber &B) {
	for (register int i = 0; i < NCIF; i++)
		if (A.C[i] != B.C[i])
			return false;
	return true;
}


int compare(BigNumber &A, BigNumber &B) {
	for (register int i = NCIF - 1; i >= 0; i--)
		if (A.C[i] != B.C[i])
			return (NCIF - i);

	return NCIF;
}
