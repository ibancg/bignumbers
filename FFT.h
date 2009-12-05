
#ifndef _FFT_H_
#define _FFT_H_

/*

 Transformadas rápidas de Fourier. (FFT e IFFT).

--------------------------------------------------------------------------------   
 (ñ) 2000 Ibán Cereijo Graña
--------------------------------------------------------------------------------
*/

#include "complex.h"

#ifndef ULI
#define ULI unsigned long int
#endif

void CreaWN();

/* ej : FFT(x, X, N) deja en X la fft de x de longitud N.
 N ha de ser potencia de 2. */
void FFT (CPX *, CPX *, ULI, ULI = 0, ULI = 0, ULI = 1);
void IFFT(CPX *, CPX *, ULI, ULI = 0, ULI = 0, ULI = 1); // trans. inversa.

#endif
