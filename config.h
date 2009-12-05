
#ifndef __CONFIG_H_
#define __CONFIG_H_

// Configuration file.

#include <math.h>

// debug mode.
#define DEBUG

// floating point precission
#define flt_t double  

#define N_DIGITS		2048  // number of digits in the format.
#define N_FRAC_DIGITS	2000  // number of fractional digits.

typedef unsigned char bcd_t; // BCD digit type.

// advanced configuration (algorithm switches).

#define FFT_MUL_ALGORITHM
#define INVERSE_DIV_ALGORITHM
#define INVERSE_NEWTON_SQRT_ALGORITHM
#define INVERSE_NEWTON_SQRT4_ALGORITHM

#endif
