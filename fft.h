#ifndef _FFT_H_
#define _FFT_H_

// Fast Fourier transforms

#include "complex.h"

void createPhaseFactors();

// computes the dft X of the signal x by a fft algorithm. The number of samples
// N is power ot 2.
void fft(Complex *x, Complex *X, unsigned long int N, unsigned long int = 0,
		unsigned long int = 0, unsigned long int = 1);

// computes the inverse fft.
void ifft(Complex *X, Complex *x, unsigned long int N, unsigned long int = 0,
		unsigned long int = 0, unsigned long int = 1);

#endif
