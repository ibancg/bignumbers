
#ifndef __CONFIG_H_
#define __CONFIG_H_

// Archivo de configuración para BIGNUM.

#include <math.h>

// modo depuración.
#define DEBUG

// precisión en coma flotante con la que trabaja el programa.
#define FLT double  

#define NCIF  1024  // número de cifras del formato.
#define NFRC  1000  // número de cifras fracionarias.

typedef unsigned char TBC; // Tipo Base Cifra.

// Configuración avanzada. (conmutadores de algoritmos).

#define ALGORITMO_MUL_FFT
#define ALGORITMO_DIV_INVERSO
#define ALGORITMO_SQRT_NEWTON_INVERSO
#define ALGORITMO_SQRT4_NEWTON_INVERSO

#endif
