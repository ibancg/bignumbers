
#include <stdio.h>

#include <stdlib.h>
#include <string.h>


#include "FFT.h"
#include "bignum.h"

int main() {

  int k;
  /*
  CreaWN();

  BigNumber DIEZ("10");
  BigNumber r;

  DIEZ.Mostrar();  
  SqrtBN(DIEZ, r);
  r.Mostrar();

  return 0;
  */
  printf("------- Cálculo del número PI (Gauss-Legendre) -------\n");
  printf("Calculando primeros iterantes...\n");
  
  CreaWN();

  BigNumber a("1");
  BigNumber b;
  BigNumber t("0.25");
  BigNumber x("1");
  BigNumber y;
  BigNumber _1p2("0.5");

  BigNumber aux1;

  SqrtBN(_1p2, b);  // b = 1/sqrt(2)

  for (k = 0; k < 20; k++) {


    //    if (ComparaBN(a, b)) break;
    printf("iteración %i\n", k);
    RestaBN(a, b, aux1);
    //    aux1.Mostrar();

    TraBN(a, y);           // y = a
    SumaBN(a, b, aux1);
    MulBN(aux1, _1p2, a);  // a = (a + b)*0.5
    MulBN(b, y, aux1);
    SqrtBN(aux1, b);       // b = sqrt(b*y)
    /*    printf("SQRT ~= ");
	  b.Mostrar();*/

    RestaBN(y, a, aux1);
    MulBN(aux1, aux1, aux1);
    MulBN(x, aux1, aux1);
    RestaBN(t, aux1, t);   // t = t - x*(y - a)^2
    SumaBN(x, x, x);       // x = 2*x;

    if (ComparaBN(a, b)) break;
  }
  
  SumaBN(t, t, t);
  SumaBN(t, t, t); // t = 4*t
  SumaBN(a, b, x);
  MulBN(x, x, x);
  DivBN(x, t, aux1); // pi = (a + b)^2/(4*t)

  printf("PI ~= ");
  aux1.Mostrar();
  printf("%u iteraciones para encontrar %u decimales de PI.\n", k, NFRC);
}
