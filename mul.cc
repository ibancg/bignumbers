
#include <string.h>

/*
--------------------------------------------------------------------------------
 (ñ) 2000 Ibán Cereijo Graña
--------------------------------------------------------------------------------

  Algoritmos de multiplicación


          AAAA.AAAA   El resultado de una multiplicación tiene el doble de
	* BBBB.BBBB  dígitos que los operandos. En un sistema de punto fijo
  -----------------  como este se perderán tantos dígitos como cifras decima-
  CCCCCCCC.CCCCCCCC  les se hayan asignado (por abajo). Por arriba se eliminan
    		     dígitos, pero no se considera pérdida de información ya
		     que si alguno de ellos fuera distinto de 0, debió haber
		     ocurrido un error de desbordamiento.
*/

#include "bignum.h"

#ifndef ALGORITMO_MUL_FFT

/* Algoritmo típico.
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
 números grandes es muy costoso.

  Si consideramos las secuencias A y B como los coeficientes de dos polino-
 mios, la convolución de estas señales representará su multiplicación. Sa-
 biendo que ningún coeficiente debe ser mayor o igual que la base del formato
 (en este caso 10), si hacemos un ajuste propagando los acarreos hacia pesos
 mayores, la secuencia resultante será la multiplicación de los dos NUMEROS
 anteriores. La multiplicación se puede hacer pues mediante una convolución
 con ajuste BCD.

  Ahora bien, un algoritmo de convolución tiene orden N^2, con lo cual no ga-
 namos nada (de hecho perder¡amos tiempo).

  Sabemos que una convolución se transforma en una multiplicación muestra por
 muestra en dominio frecuencial. Así, podemos transformar las señales (DFT),
 multiplicarlas en frecuencia muestra a muestra y aplicar una transformación
 inversa.

  Calcular una DFT (Discrete Fourier Tranform) se haría normalmente mediante
 un algoritmo de orden N^2, pero aplicando un algoritmo FFT (Fast Fourier
 Transform), el orden del algoritmo se reduce a N*log2(N), con lo que resulta
 mucho mejor para N grande.

  Optimizando:

  La FFT de una señal será compleja aunque no lo sea la señal de entrada, y
 sabemos que las dos señales de entrada (números a multiplicar) son señales
 reales. Así, podemos construir una sola señal que contenga la información de
 las dos señales en su parte real y en su parte imaginaria.

    x[n] = x1 + i*x2[n] --> Calculamos su FFT --> X[n]

  Podemos extraer, usando las propiedades de la FFT:

	X1[n] =     1/2*(X[n] + conj(X[-n mod N]))
	X2[n] = 1/(2*i)*(X[n] - conj(X[-n mod N]))

  De esta forma sólo es necesario hacer una FFT y no 2.

------------------------------------------------------------------------------
*/

#include "complex.h"
#include "FFT.h"

// Implementación del algoritmo. Todos solapables. p.ej : MulBN(X, X, X)
void mul(BigNumber &A, BigNumber &B, BigNumber &C)
{
  static CPX         BC1[NCIF << 1];
  static CPX         BC2[NCIF << 1];
  register int       i;
  CPX                Xi, Xmi, X1, X2, X3;


  // paso 1: construyo una señal compleja que contenga información de las 2.
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

  // paso 4: destransformo la señal.
  IFFT(BC1, BC2, (NCIF << 1));

  unsigned long int  c, ci;
  FLT                x;

  // paso 5: limpio y ajusto la señal.
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
