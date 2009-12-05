
#include <string.h>

/*
--------------------------------------------------------------------------------
 (�) 2000 Ib�n Cereijo Gra�a
--------------------------------------------------------------------------------

  Algoritmos de multiplicaci�n


          AAAA.AAAA   El resultado de una multiplicaci�n tiene el doble de
	* BBBB.BBBB  d�gitos que los operandos. En un sistema de punto fijo
  -----------------  como este se perder�n tantos d�gitos como cifras decima-
  CCCCCCCC.CCCCCCCC  les se hayan asignado (por abajo). Por arriba se eliminan
    		     d�gitos, pero no se considera p�rdida de informaci�n ya
		     que si alguno de ellos fuera distinto de 0, debi� haber
		     ocurrido un error de desbordamiento.
*/

#include "bignum.h"

#ifndef ALGORITMO_MUL_FFT

/* Algoritmo t�pico.
 Todos se pueden solapar. ej : MulBN(X, X, X) */
void mul(BigNumber &A, BigNumber &B, BigNumber &C)
{
  int               i, j;
  unsigned long int r, c;
  static TBC        R[NCIF + NFRC];


  memset(R, 0, (NCIF + NFRC)*sizeof(TBC));
  
  for (i = 0; i < NCIF; i++) {
    for (j = 0, c = 0; (j < NCIF) && ((j + i) < (NCIF + NFRC)); j++) {
      
      r = R[j + i] + c + A.C[j]*B.C[i];
      c = (r/10);
      R[j + i] = (r - 10*c);
    }

    for (; (j + i) < (NCIF + NFRC); j++) { // ajusto lo que me queda.
      r = R[j + i] + c;
      c = (r/10);
      R[j + i] = (r - 10*c);
    }
  }
  
  C.isPositive = !(A.isPositive ^ B.isPositive);
  memcpy(C.C, &R[NFRC], NCIF*sizeof(TBC));
}


#else

/*
------------------------------------------------------------------------------

  El problema del algoritmo arriba propuesto es que tiene orden N^2. Para
 n�meros grandes es muy costoso.

  Si consideramos las secuencias A y B como los coeficientes de dos polino-
 mios, la convoluci�n de estas se�ales representar� su multiplicaci�n. Sa-
 biendo que ning�n coeficiente debe ser mayor o igual que la base del formato
 (en este caso 10), si hacemos un ajuste propagando los acarreos hacia pesos
 mayores, la secuencia resultante ser� la multiplicaci�n de los dos NUMEROS
 anteriores. La multiplicaci�n se puede hacer pues mediante una convoluci�n
 con ajuste BCD.

  Ahora bien, un algoritmo de convoluci�n tiene orden N^2, con lo cual no ga-
 namos nada (de hecho perder�amos tiempo).

  Sabemos que una convoluci�n se transforma en una multiplicaci�n muestra por
 muestra en dominio frecuencial. As�, podemos transformar las se�ales (DFT),
 multiplicarlas en frecuencia muestra a muestra y aplicar una transformaci�n
 inversa.

  Calcular una DFT (Discrete Fourier Tranform) se har�a normalmente mediante
 un algoritmo de orden N^2, pero aplicando un algoritmo FFT (Fast Fourier
 Transform), el orden del algoritmo se reduce a N*log2(N), con lo que resulta
 mucho mejor para N grande.

  Optimizando:

  La FFT de una se�al ser� compleja aunque no lo sea la se�al de entrada, y
 sabemos que las dos se�ales de entrada (n�meros a multiplicar) son se�ales
 reales. As�, podemos construir una sola se�al que contenga la informaci�n de
 las dos se�ales en su parte real y en su parte imaginaria.

    x[n] = x1 + i*x2[n] --> Calculamos su FFT --> X[n]

  Podemos extraer, usando las propiedades de la FFT:

	X1[n] =     1/2*(X[n] + conj(X[-n mod N]))
	X2[n] = 1/(2*i)*(X[n] - conj(X[-n mod N]))

  De esta forma s�lo es necesario hacer una FFT y no 2.

------------------------------------------------------------------------------
*/

#include "complex.h"
#include "FFT.h"

// Implementaci�n del algoritmo. Todos solapables. p.ej : MulBN(X, X, X)
void mul(BigNumber &A, BigNumber &B, BigNumber &C)
{
  static CPX         BC1[NCIF << 1];
  static CPX         BC2[NCIF << 1];
  register int       i;
  CPX                Xi, Xmi, X1, X2, X3;


  // paso 1: construyo una se�al compleja que contenga informaci�n de las 2.
  for (i = 0; i < NCIF; i++) {
    BC1[i].r = A.C[i]; // parte real
    BC1[i].i = B.C[i]; // parte imaginaria.
  }
  memset(&BC1[NCIF], 0, NCIF*sizeof(CPX)); // limpio la parte alta.

  // paso 2: la transformo.
  FFT(BC1, BC2, (NCIF << 1));
  
  // paso 3: multiplico en frecuencia.
  for (i = 0; i < (NCIF << 1); i++) {

    // hay que extraer las transformadas individuales de la trans. conjunta.
    Xi  = BC2[i];
    Xmi = BC2[(-i) & ((NCIF << 1) - 1)];
    Xmi.i = -Xmi.i; // conjugado

    X1.r = Xi.r + Xmi.r;
    X1.i = Xi.i + Xmi.i; // X1 = Xi + Xmi
    X2.r = Xi.r - Xmi.r;
    X2.i = Xi.i - Xmi.i; // X2 = Xi - Xmi

    // ahora multiplico muestra por muestra.
    X3.r = X1.r*X2.r - X1.i*X2.i;
    X3.i = X1.r*X2.i + X1.i*X2.r; // X3 = X1*X2;

    BC1[i].r =  0.25*X3.i;
    BC1[i].i = -0.25*X3.r;
  }

  // paso 4: destransformo la se�al.
  IFFT(BC1, BC2, (NCIF << 1));

  unsigned long int  c, ci;
  FLT                x;

  // paso 5: limpio y ajusto la se�al.
  for (i = 0, c = 0; i < NCIF + NFRC; i++) {

    x = BC2[i].r; // quito la parte imaginaria.

    if ((x - floor(x)) < 0.5) ci = (unsigned long int)(c + floor(x));
    else ci = (unsigned long int)(c + ceil(x)); // redondeo.

    c = (ci/10);                                // propago el acarreo
    if (i >= NFRC) C.C[i - NFRC] = (ci - 10*c); // ci % 10
  }
   
  C.isPositive = !(A.isPositive ^ B.isPositive);   
}

#endif
