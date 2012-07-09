
#include <stdio.h>
#include <stdlib.h>

#include "bignum.hh"

void main() {

  // Comprobación de una terna de Fermat.

  BigNumber X1, X2, X3, AX1, AX2, AX3, AX;
   
  X1 = BigNumber("1841");
  X2 = BigNumber("1782");
  X3 = BigNumber("1922");
   
  printf("Calculando...\n");
   
  AX1 = X1;
  AX2 = X2;
  AX3 = X3;
 
  for (int i = 1; i < 12; i++) {

    mul(AX1, X1, AX);
    AX1 = X1;
    mul(AX2, X2, AX);
    AX2 = X2;
    mul(AX3, X3, AX);
    AX3 = X3;
  }   
   
  add(AX1, AX2, AX);
  AX.show();
  AX3.show();
}
