
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bignum.h"

void Inverso(BigNumber &, BigNumber &);
void SqrtmBN(BigNumber &, BigNumber &);

void main() {

  char     s[80];
  double   d;

  BigNumber X, Y, Z;
  int       i, j;
  
  //  FLT2BN(1.1e-20, X);
  X = BigNumber("7");
  Y = BigNumber();
    
  Inverso(X, Z);
  Z.show();
  /*		DivBN(X, Y, Z);
		Z.Mostrar();*/

  /*	X.Mostrar();
	Y.Mostrar();
	SqrtBN(X, Y);
	Y.Mostrar();*/
}

//-------------------------------------------------------------------

void Inverso(BigNumber &A, BigNumber &x)
{
  static BigNumber x1, x2;
  static BigNumber DOS("2");
  int              i;
  flt_t              d;
  int              ipc = 0; // índice primera cifra.
  bool             stop;

  x.isPositive = A.isPositive;
  memset(x.digits, 0, BigNumber::N_DIGITS*sizeof(bcd_t)); // limpiamos B.

  // cuento el número de cifras enteras.
  /*  for (i = NCIF - 1; (i >= 0) && !ipc; i--) if (A.C[i]) ipc = i;
      x.C[NFRC - (ipc - NFRC) - 1] = 1;*/

  d = bigNumber2Flt(A);
  flt2BigNumber(1.0/d, x);

  for (int k = 0;; k++) {
    
    mul(A, x, x1);
    sub(DOS, x1, x2);
    mul(x, x2, x1);

    for (stop = true, i = 0; (i < BigNumber::N_DIGITS) && stop; i++) 
      if (x.digits[i] != x1.digits[i]) stop = false;

    if (stop) break;

    copy(x1, x);

    printf("INVERSO : iteración %i\n", k);
    x.show();
  }

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
    
    printf("SQRT : iteracion = %i\n", k);
    x.show();
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
