
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
   
  TraBN(X1, AX1);
  TraBN(X2, AX2);
  TraBN(X3, AX3);
 
  for (int i = 1; i < 12; i++) {

    MulBN(AX1, X1, AX);
    TraBN(AX, AX1);
    MulBN(AX2, X2, AX);
    TraBN(AX, AX2);
    MulBN(AX3, X3, AX);
    TraBN(AX, AX3);
  }   
   
  SumaBN(AX1, AX2, AX);
  AX.Mostrar();
  AX3.Mostrar();
}
