
#include <stdio.h>
//#include <conio.h>
#include <stdlib.h>
#include <string.h>

#include "bignum.h"

void SqrtmBN(BigNumber &A, BigNumber &x);

void main() {

  
  printf("------- Cálculo del número PI -------\n");
  printf("Calculando primeros iterantes...\n");
  
  BigNumber UNO("1");
  BigNumber ACP2("2"); // acumulador de potencias de 2.
  BigNumber y, a, x1, x2, x3, pi, pio;
  bool      stop, z;
  int       i, j;

  pi = BigNumber("0");
  
  SqrtmBN(ACP2, x2);
  sub(x2, UNO, y); // y0 = sqrt(2) - 1
  
  add(x2, x2, x2); // 2*sqrt(2)
  add(x2, x2, x2); // 4*sqrt(2)
  
  a = BigNumber("6");
  sub(a, x2, a);  // a0 = 6 - 4*sqrt(2)
  printf("...OK\n");
  
  for (i = 0;; i++) {
    
    div(UNO, a, pio);
    //    gotoxy(1, 4);

    // paramos cuando dos iterantes consecutivos coinciden.
    for (stop = true, j = BigNumber::N_DIGITS - 1; stop && (j >= 0); j--)
      if (pio.digits[j] != pi.digits[j]) stop = false;

    printf("iteración %uª : %u decimales encontrados\n", i + 1, BigNumber::N_FRAC_DIGITS - j - 1);

    if (stop) break;
    copy(pio, pi);

    mul(y, y, x1);   // y^2
    mul(x1, x1, x2); // y^4
    sub(UNO, x2, x1); // 1 - y^4
	
    /*    Sqrt4BN(x1, x2); // (1 - y^4)^(1/4);
	  RestaBN(UNO, x2, x1); // (1 - (1 - y^4)^(1/4))
	  SumaBN(UNO, x2, x3);  // (1 + (1 - y^4)^(1/4))
	  DivBN(x1, x3, y);     // (1 - (1 - y^4)^(1/4))/(1 + (1 - y^4)^(1/4))*/

    SqrtmBN(x1, x2); // (1 - y^4)^(1/2);
    SqrtmBN(x2, x1); // (1 - y^4)^(1/4);
    sub(UNO, x1, x2); // (1 - (1 - y^4)^(1/4))
    add(UNO, x1, x3);  // (1 + (1 - y^4)^(1/4))
    div(x2, x3, y);     // (1 - (1 - y^4)^(1/4))/(1 + (1 - y^4)^(1/4))

    add(y, UNO, x1); // (1 + y)
    mul(x1, x1, x2);  // (1 + y)^2
    mul(x2, x2, x3);  // (1 + y)^4
    mul(x3, a, x2);   // (1 + y)^4*a

    add(ACP2, ACP2, ACP2); // 2*ACP2
    add(ACP2, ACP2, ACP2); // 4*ACP2

    mul(y, y, x3); // y^2
    add(x1, x3, x1); // (1 + y + y^2)
    mul(y, x1, x3);     // y*(1 + y + y^2)	
    mul(ACP2, x3, x1); // 2^(2*i + 1)*y*(1 + y + y^2)
    sub(x2, x1, a); // (1 + y)^4*a - 2^(2*i + 1)*y*(1 + y + y^2)
    
  }
  
  printf("PI ~= ");
  pi.show();
  printf("%u iteraciones para encontrar %u cifras decimales de PI.\n", i, BigNumber::N_FRAC_DIGITS);
      
}

// Cálculo de la raíz cuadrada por el método de Newton. (F(x) = x^2 - 1/A).
// A y x NO solapables.
void SqrtmBN(BigNumber &A, BigNumber &x)
{
  static BigNumber    x2, Fx, DFx, xo;
  static BigNumber    _1p2("0.5");
  static BigNumber    TRES("3");
  int                 i, n;
  bool                stop;

   
#ifdef DEBUG
  int _x, _y;
  _x = wherex();
  _y = wherey();
  textcolor(14);
  gotoxy(__X, 1);
  cprintf(" sqrt2 ");
  __X += 10;
#endif

  if (!A.isPositive) {
    printf("ERROR: raíz compleja\n");
    exit(255);
  }

  /* Si el orden de magnitud del número es n, empiezo a iterar en 10^(n/2)*/
  for (i = BigNumber::N_DIGITS - 1, n = -1; (i >= 0) && (n == -1); i--)
    if (A.digits[i]) n = i;

  xo.isPositive = true;
  memset(xo.digits, 0, BigNumber::N_DIGITS*sizeof(bcd_t));

  if (n == -1) { // es un 0
    copy(xo, x);
    return;
  }

  n = BigNumber::N_FRAC_DIGITS - (n - BigNumber::N_FRAC_DIGITS + 1)/2;
  xo.digits[n] = 1;
     
  //xo.Mostrar();
  //  getchar();

  for (int k = 0;; k++) {

    copy(xo, x);
    mul(x, x, x2);
    mul(x2, A, Fx);              // A*x^2
    sub(Fx, TRES, DFx, false); // (A*x^2 - 3)
    mul(_1p2, DFx, Fx);          // 0.5*(A*x^2 - 3)
    mul(x, Fx, xo);              // 0.5*x*(A*x^2 - 3)
    xo.isPositive = !xo.isPositive;

    // Si se repite el iterante, paramos.
    for (stop = true, i = 0; (i < BigNumber::N_DIGITS) && stop; i++) 
      if (x.digits[i] != xo.digits[i]) stop = false;
    
    //    printf("SQRT : iteracion = %i\n", k);
    if (stop) break;
  }
  
  mul(A, xo, x);

#ifdef DEBUG
  gotoxy(__X, 1);
  cprintf("       ");
  __X -= 10;
  textcolor(7);
  gotoxy(_x, _y);
#endif

  // Por si el método de Newton me converge a la solución negativa
  x.isPositive = true;
  
  // Si F(x) <= 0 -> x^2 <= A, podemos salir. si no, hay que restar 1.
  if (!Fx.isPositive) return; // F(x) <= 0. (signo forzado en 0).

  char c = 1;
  
  // ajustamos (debido a esta resta) hasta donde tengamos que hacerlo.
  for (i = 0;; i++) {
    x.digits[i] -= c;
    c = (((char)x.digits[i]) < 0) ? 1 : 0;
    if (!c) break; // si no hay acarreo salgo.
    x.digits[i] += 10*c;
  }
}
