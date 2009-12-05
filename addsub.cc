
#include <string.h>

#include "bignum.h"

/* Suma. Todos se pueden solapar. ej: SumaBN(x, x, x);
 en caso de que el resultado sea 0, el signo se fuerza con piz */
void SumaBN(BigNumber &A, BigNumber &B, BigNumber &C, bool piz)
{
  int  i;
  char r;
  char c = 0; // acarreo.

  if (!(A.positivo ^ B.positivo)) { // mismo signo.
	   
    for (i = 0; i < NCIF; i++) {

      r = c + A.C[i] + B.C[i];
      c = (r > 9) ? 1 : 0;
      C.C[i] = (r - 10*c);  // r % 10
    }
	 
    C.positivo = A.positivo; // = B.positivo
    return;
  }
   
  // Distinto signo.
  BigNumber *M, *m; // número de mayor y menor módulo.

  M = NULL;
  
  for (i = NCIF - 1; i >= 0; i--) {
	   
    if (A.C[i] == B.C[i]) continue;
	   
    if (A.C[i] > B.C[i]) {
      M = &A;
      m = &B;
    } else {
      M = &B;
      m = &A;   
    }
    break;
  }

  if (!M) { // son iguales.	  
    memset(C.C, 0, NCIF*sizeof(TBC));
    C.positivo = piz;
    return;
  }  	  
  
  // Al de mayor módulo le resto el de menor.   
  for (i = 0; i < NCIF; i++) {

    r = M->C[i] - (m->C[i] + c);
    c = (r < 0) ? 1 : 0;
    C.C[i] = (r + 10*c);
  }

  // Si el positivo es el de mayor módulo el resultado tiene signo positivo.
  C.positivo = ((A.positivo) && (M == &A)) || ((B.positivo) && (M == &B));
}

//------------------------------------------------------------------------

// Resta. Todos se pueden solapar.
void RestaBN(BigNumber &A, BigNumber &B, BigNumber &C, bool piz)
{
  int  i;
  char r;
  char c = 0; // acarreo.

  if (!(A.positivo ^ !B.positivo)) { // mismo signo.
	   
    for (i = 0; i < NCIF; i++) {

      r = c + A.C[i] + B.C[i];
      c = (r > 9) ? 1 : 0;
      C.C[i] = (r - 10*c);  // r % 10
    }
	 
    C.positivo = A.positivo; // = B.positivo
    return;
  }
   
  // Distinto signo.
  BigNumber *M, *m; // número de mayor y menor módulo.

  M = NULL;
  
  for (i = NCIF - 1; i >= 0; i--) {
	   
    if (A.C[i] == B.C[i]) continue;
	   
    if (A.C[i] > B.C[i]) {
      M = &A;
      m = &B;
    } else {
      M = &B;
      m = &A;   
    }
    break;
  }

  if (!M) { // son iguales.
	  
    memset(C.C, 0, NCIF*sizeof(TBC));
    C.positivo = piz;
    return;
  }  	  
  
  // Al de mayor módulo le resto el de menor.   
  for (i = 0; i < NCIF; i++) {

    r = M->C[i] - (m->C[i] + c);
    c = (r < 0) ? 1 : 0;
    C.C[i] = (r + 10*c);
  }

  // Si el positivo es el de mayor módulo el resultado tiene signo positivo.
  C.positivo = ((A.positivo) && (M == &A)) || ((!B.positivo) && (M == &B));
}
