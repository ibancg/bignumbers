/*
 BigNumbers - Arbitrary precision arithmetic
 Copyright 2000-2010, Ibán Cereijo Graña <ibancg at gmail dot com>

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

#include "config.h"
#include "complex.h"
#include "fft.h"

// Phase factors.
Complex* WN;

// Optimization: computes the phase factor table WN[i] = exp(-j*k), with
// k = 0..pi (N samples)
void createPhaseFactors() {

	WN = new Complex[N_DIGITS];
	flt_t alpha;

	for (int i = 0; i < N_DIGITS; i++) {
		alpha = -i * M_PI / N_DIGITS;
		WN[i].r = cos(alpha);
		WN[i].i = sin(alpha);
	}
}

void destroyPhaseFactors() {
	delete[] WN;
}

void fft(Complex *x, Complex *X, unsigned long int N, unsigned long int offset,
		unsigned long int d1, unsigned long int step) {
	Complex X1, X2;
	unsigned long int Np2 = (N >> 1); // N/2
	register unsigned long int a, b, c, q;

	if (N == 2) { // Butterfly for N = 2;

		X1 = x[offset];
		X2 = x[offset + step];
		X[d1].r = X1.r + X2.r;
		X[d1].i = X1.i + X2.i; // X[d1] = X1 + X2
		X[d1 + Np2].r = X1.r - X2.r;
		X[d1 + Np2].i = X1.i - X2.i; // X[d1 + Np2] = X1 - X2
		return;
	}

	fft(x, X, Np2, offset, d1, step << 1);
	fft(x, X, Np2, offset + step, d1 + Np2, step << 1);

	for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

		a = q + d1;
		b = a + Np2;

		X1 = X[a];
		X2.r = X[b].r * WN[c].r - X[b].i * WN[c].i;
		X2.i = X[b].r * WN[c].i + X[b].i * WN[c].r; // X2 = X[b]*WN[c]

		X[a].r = X1.r + X2.r;
		X[a].i = X1.i + X2.i; // X[a] = X1 + X2
		X[b].r = X1.r - X2.r;
		X[b].i = X1.i - X2.i; // X[b] = X1 - X2
	}
}

void ifft(Complex *X, Complex *x, unsigned long int N,
		unsigned long int offset, unsigned long int d1, unsigned long int step) {
	Complex x1, x2;
	unsigned long int Np2 = (N >> 1); // N/2
	flt_t _1pN = 1.0 / N;
	register unsigned long int a, b, c, q;

	if (N == 2) { // Butterfly for N = 2;

		x1 = X[offset];
		x2 = X[offset + step];
		x[d1].r = x1.r + x2.r;
		x[d1].i = x1.i + x2.i; // x[d1] = x1 + x2
		x[d1 + Np2].r = x1.r - x2.r;
		x[d1 + Np2].i = x1.i - x2.i; // x[d1 + Np2] = x1 - x2

		return;
	}

	ifft(X, x, Np2, offset, d1, step << 1);
	ifft(X, x, Np2, offset + step, d1 + Np2, step << 1);

	for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

		a = q + d1;
		b = a + Np2;

		x1 = x[a];
		x2.r = x[b].r * WN[c].r + x[b].i * WN[c].i;
		x2.i = x[b].i * WN[c].r - x[b].r * WN[c].i; // x2 = x[b]*WN*[c]

		x[a].r = x1.r + x2.r;
		x[a].i = x1.i + x2.i; // x[a] = x1 + x2
		x[b].r = x1.r - x2.r;
		x[b].i = x1.i - x2.i; // x[b] = x1 - x2
	}

	if (step != 1)
		return;

	_1pN = 1.0 / N;
	for (q = 0; q < N; q++) {
		x[q].r = x[q].r * _1pN;
		x[q].i = x[q].i * _1pN; // x[q] = x[q]/N
	}
}
