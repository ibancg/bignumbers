/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2010, Ib�n Cereijo Gra�a <ibancg at gmail dot com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _FFT_H_
#define _FFT_H_

// Fast Fourier transforms

#include <complex>
#include <vector>

// initializes the libreary
void createTwiddleFactors(unsigned long int N);

// computes the dft X of the signal x by a fft algorithm. The number of samples
// N must be a power ot 2.
void fft(const std::vector<std::complex<double> >& x,
		std::vector<std::complex<double> >& X, unsigned long int N,
		unsigned long int = 0, unsigned long int = 0, unsigned long int = 1);

// computes the inverse fft.
void ifft(const std::vector<std::complex<double> >& X,
		std::vector<std::complex<double> >& x, unsigned long int N,
		unsigned long int = 0, unsigned long int = 0, unsigned long int = 1);

#endif
