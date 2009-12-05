
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

  sqrt(_1p2, b);  // b = 1/sqrt(2)

  for (k = 0; k < 20; k++) {


    //    if (ComparaBN(a, b)) break;
    printf("iteración %i\n", k);
    sub(a, b, aux1);
    //    aux1.Mostrar();

    copy(a, y);           // y = a
    add(a, b, aux1);
    mul(aux1, _1p2, a);  // a = (a + b)*0.5
    mul(b, y, aux1);
    sqrt(aux1, b);       // b = sqrt(b*y)
    /*    printf("SQRT ~= ");
	  b.Mostrar();*/

    sub(y, a, aux1);
    mul(aux1, aux1, aux1);
    mul(x, aux1, aux1);
    sub(t, aux1, t);   // t = t - x*(y - a)^2
    add(x, x, x);       // x = 2*x;

    if (equals(a, b)) break;
  }
  
  add(t, t, t);
  add(t, t, t); // t = 4*t
  add(a, b, x);
  mul(x, x, x);
  div(x, t, aux1); // pi = (a + b)^2/(4*t)

  printf("PI ~= ");
  aux1.show();
  printf("%u iteraciones para encontrar %u decimales de PI.\n", k, NFRC);
}
