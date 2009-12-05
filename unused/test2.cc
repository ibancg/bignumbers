
#include <stdio.h>
#include <stdlib.h>

#include "bignum.hh"

void main() {

	BigNumber X, Y, Z;
	int       i, j;
	
	X = BigNumber("8798476890740666256565");
	Y = BigNumber("7852762856720958672665");
	
/*      	for (i = 0; i < NCIF; i++) {
	  X.C[i] = (unsigned char)(10.0*random()/RAND_MAX);
	  Y.C[i] = (unsigned char)(10.0*random()/RAND_MAX);
	}
	
	X.positivo = true;
	Y.positivo = true;*/

	X.Mostrar();	
	Y.Mostrar();	
	printf("Multiplicando...\n");
	for (i = 0; i < 30; i++) MulBN(X, Y, Z);
	Z.Mostrar();	
/*	printf("Dividiendo...\n");
	for (i = 0; i < 100; i++) DivBN(Z, Y, X);*/
	//	X.Mostrar();
}
