
#ifndef __CONFIG_H_
#define __CONFIG_H_

// Archivo de configuraci�n para BIGNUM.

#include <math.h>

// modo depuraci�n.
#define DEBUG

// precisi�n en coma flotante con la que trabaja el programa.
#define FLT double  

#define NCIF  1024  // n�mero de cifras del formato.
#define NFRC  1000  // n�mero de cifras fracionarias.

typedef unsigned char TBC; // Tipo Base Cifra.

// Configuraci�n avanzada. (conmutadores de algoritmos).

#define ALGORITMO_MUL_FFT
#define ALGORITMO_DIV_INVERSO
#define ALGORITMO_SQRT_NEWTON_INVERSO
#define ALGORITMO_SQRT4_NEWTON_INVERSO

#endif
