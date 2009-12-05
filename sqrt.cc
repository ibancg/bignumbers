
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "bignum.h"

using namespace std;

#ifdef ALGORITMO_SQRT_NEWTON_INVERSO

// Cálculo de la raíz cuadrada por el método de Newton. (F(x) = x^2 - 1/A).
// A y x NO solapables.
void SqrtBN(BigNumber &A, BigNumber &x)
{
  static BigNumber    xo;
  static BigNumber    _1p2("0.5");
  static BigNumber    TRES("3");

# ifdef DEBUG
  cout << "SQRT";
  cout.flush();
# endif

  if (!A.positivo) {
    printf("ERROR: raíz compleja.\n");
    exit(255);
  }

  xo.positivo = true;
  memset(xo.C, 0, NCIF*sizeof(TBC));
  
  // Si el orden de magnitud del número es n, empiezo a iterar en 10^(-n/2)
  xo.C[NFRC - (IndicePC(A) - NFRC + 1)/2] = 1;
  /* FLT d
     BN2FLT(A, d);
     FLT2BN(sqrt(d), xo);*/
  
  for (int k = 0;; k++) {
    
    TraBN(xo, x);
    MulBN(x, x, xo);
    MulBN(xo, A, xo);             // A*x^2
    RestaBN(xo, TRES, xo, false); // (A*x^2 - 3)
    MulBN(_1p2, xo, xo);          // 0.5*(A*x^2 - 3)
    MulBN(x, xo, xo);             // 0.5*x*(A*x^2 - 3)
    xo.positivo = !xo.positivo;

#   ifdef DEBUG
    cout << '.';
    cout.flush();
#   endif
    
    if (Compara2BN(x, xo) >= (NCIF - 2)) break; // Si se repite el iterante, paramos.
  }
  
# ifdef DEBUG
  cout << endl;
# endif

  xo.positivo = true; // por si converge a la solución negativa.
  MulBN(A, xo, x); // el método ha convergido al inverso.
}

#else
//------------------------------------------------------------------------

// Cálculo de la raíz cuadrada por el método de Newton. (F(x) = x^2 - A).
// A y x NO solapables.
void SqrtBN(BigNumber &A, BigNumber &x)
{
  static BigNumber    x2, Fx, DFx, xo;

# ifdef DEBUG
  cout << "SQRT";
  cout.flush();
# endif

  if (!A.positivo) {
    printf("ERROR: raíz compleja.\n");
    exit(255);
  }

  xo.positivo = true;
  memset(xo.C, 0, NCIF*sizeof(TBC));
  
  /* Si el orden de magnitud del número es n, empiezo a iterar en 10^(n/2)*/
  xo.C[NFRC + (IndicePC(A) - NFRC)/2] = 1;
  /*  FLT d;
      BN2FLT(A, d);
      FLT2BN(pow(d, 0.5), xo); */
  
  for (;;) {
	
    TraBN(xo, x);
    MulBN(x, x, x2);
    RestaBN(x2, A, Fx, false); // F(x) = x^2 - A
	 
    SumaBN(x, x, DFx);  // F'(x) = 2*x
    DivBN(Fx, DFx, x2); // F(x)/F'(x) 
    RestaBN(x, x2, xo); // x - F(x)/F'(x) 

#   ifdef DEBUG
    cout << '.';
    cout.flush();
#   endif
    
    if (ComparaBN(x, xo)) break; // Si se repite el iterante, paramos.
  }
  
# ifdef DEBUG
  cout << endl;
# endif

  // Por si el método de Newton me converge a la solución negativa
  x.positivo = true;
  
  // Si F(x) <= 0 -> x^2 <= A, podemos salir. si no, hay que restar 1.
  if (!Fx.positivo) return; // F(x) <= 0. (signo forzado en 0).

  char c = 1;
  
  // ajustamos (debido a esta resta) hasta donde tengamos que hacerlo.
  for (int i = 0;; i++) {
    x.C[i] -= c;
    c = (((char)x.C[i]) < 0) ? 1 : 0;
    if (!c) break; // si no hay acarreo salgo.
    x.C[i] += 10*c;
  }
}

#endif

//---------------------------- RAIZ CUARTA ----------------------------

#ifdef ALGORITMO_SQRT4_NEWTON_INVERSO

// Cálculo de la raíz cuarta por el método de Newton. (F(x) = x^4 - 1/A).
// A y x NO solapables.
void Sqrt4BN(BigNumber &A, BigNumber &x)
{
  static BigNumber    xo;
  static BigNumber    _1p4("0.25");
  static BigNumber    CINCO("5");

  if (!A.positivo) {
    printf("ERROR: raíz compleja.\n");
    exit(255);
  }

# ifdef DEBUG
  cout << "SQRT4";
  cout.flush();
# endif

  xo.positivo = true;
  memset(xo.C, 0, NCIF*sizeof(TBC));
      
  /* Si el orden de magnitud del número es n, empiezo a iterar en 10^(n/4)*/
  xo.C[NFRC - (IndicePC(A) - NFRC + 1)/4] = 1;
  /* FLT d;
     BN2FLT(A, d);
     FLT2BN(pow(d, 0.25), xo);*/

  for (int k = 0;; k++) {
    
    TraBN(xo, x);
    MulBN(x, x, xo);
    MulBN(xo, xo, xo);             // x^4
    MulBN(xo, A, xo);              // A*x^4
    RestaBN(xo, CINCO, xo, false); // (A*x^2 - 5)
    MulBN(_1p4, xo, xo);           // 0.25*(A*x^4 - 5)
    MulBN(x, xo, xo);              // 0.25*x*(A*x^4 - 5)
    xo.positivo = !xo.positivo;

#   ifdef DEBUG
    cout << '.';
    cout.flush();
#   endif
    
    if (Compara2BN(x, xo) >= (NCIF - 2)) break; // Si se repite el iterante, paramos.
  }
  
# ifdef DEBUG
  cout << endl;
# endif

  xo.positivo = true; // por si converge a la solución negativa.
  InvBN(xo, x);
}

#else 
//------------------------------------------------------------------------

// Cálculo de la raíz cuarta por el método de Newton. (F(x) = x^4 - A).
// A y x NO solapables.
void Sqrt4BN(BigNumber &A, BigNumber &x)
{
  static BigNumber    x1, x2, Fx, DFx, xo;
 
  if (!A.positivo) {
    printf("ERROR: raíz compleja.\n");
    exit(255);
  }

  xo.positivo = true;
  memset(xo.C, 0, NCIF*sizeof(TBC));

  /* Si el orden de magnitud del n£mero es n, empiezo a iterar en 10^(n/4)*/
  xo.C[NFRC + (IndicePC(A) - NFRC)/4] = 1;
  /*  FLT d;
      BN2FLT(A, d);
      FLT2BN(pow(d, 0.25), xo); */

  for (;;) {

    TraBN(xo, x);
    MulBN(x, x, x1);  // x^2
    MulBN(x, x1, x1); // x^3
    MulBN(x, x1, x2); // x^4
    RestaBN(x2, A, Fx, false); // F(x) = x^4 - A
	 
    SumaBN(x1, x1, x1);  // 2*x^3
    SumaBN(x1, x1, DFx); // F'(x) = 4*x^3
    DivBN(Fx, DFx, x2); // F(x)/F'(x) 
    RestaBN(x, x2, xo); // x - F(x)/F'(x) 
	
    if (ComparaBN(x, xo)) break; // Si se repite el iterante, paramos.
  }
  
  // Por si el método de Newton me converge a la solución negativa
  x.positivo = true;
  
  // Si F(x) <= 0 -> x^2 <= A, podemos salir. si no, hay que restar 1.
  if (!Fx.positivo) return; // F(x) <= 0. (signo forzado en 0).

  char c = 1;
  
  // ajustamos (debido a esta resta) hasta donde tengamos que hacerlo.
  for (int i = 0;; i++) {
    x.C[i] -= c;
    c = (((char)x.C[i]) < 0) ? 1 : 0;
    if (!c) break; // si no hay acarreo salgo.
    x.C[i] += 10*c;
  }
}

#endif
