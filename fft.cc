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
#include "fft.h"
#include "bignum.h"

// Phase factors.
std::complex<double>* WN;

// Optimization: computes the phase factor table WN[i] = exp(-j*k), with
// k = 0..pi (N samples)
void createPhaseFactors() {

	WN = new std::complex<double>[BigNumber::N_DIGITS];
	double alpha;

	for (int i = 0; i < BigNumber::N_DIGITS; i++) {
		alpha = -i * M_PI / BigNumber::N_DIGITS;
		WN[i] = std::complex<double>(cos(alpha), sin(alpha));
	}
}

void destroyPhaseFactors() {
	delete[] WN;
}

void fft(const std::vector<std::complex<double> >& x,
		std::vector<std::complex<double> >& X, unsigned long int N,
		unsigned long int offset, unsigned long int d1,
		unsigned long int step) {
	std::complex<double> X1, X2;
	unsigned long int Np2 = (N >> 1); // N/2
	register unsigned long int a, b, c, q;

	if (N == 2) { // Butterfly for N = 2;

		X1 = x[offset];
		X2 = x[offset + step];
		X[d1].real() = X1.real() + X2.real();
		X[d1].imag() = X1.imag() + X2.imag(); // X[d1] = X1 + X2
		X[d1 + Np2].real() = X1.real() - X2.real();
		X[d1 + Np2].imag() = X1.imag() - X2.imag(); // X[d1 + Np2] = X1 - X2
		return;
	}

	fft(x, X, Np2, offset, d1, step << 1);
	fft(x, X, Np2, offset + step, d1 + Np2, step << 1);

	for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

		a = q + d1;
		b = a + Np2;

		X1 = X[a];
		X2.real() = X[b].real() * WN[c].real() - X[b].imag() * WN[c].imag();
		X2.imag() = X[b].real() * WN[c].imag() + X[b].imag() * WN[c].real(); // X2 = X[b]*WN[c]

		X[a].real() = X1.real() + X2.real();
		X[a].imag() = X1.imag() + X2.imag(); // X[a] = X1 + X2
		X[b].real() = X1.real() - X2.real();
		X[b].imag() = X1.imag() - X2.imag(); // X[b] = X1 - X2
	}
}

void ifft(const std::vector<std::complex<double> >& X,
		std::vector<std::complex<double> >& x, unsigned long int N,
		unsigned long int offset, unsigned long int d1,
		unsigned long int step) {
	std::complex<double> x1, x2;
	unsigned long int Np2 = (N >> 1); // N/2
	double _1pN = 1.0 / N;
	register unsigned long int a, b, c, q;

	if (N == 2) { // Butterfly for N = 2;

		x1 = X[offset];
		x2 = X[offset + step];
		x[d1].real() = x1.real() + x2.real();
		x[d1].imag() = x1.imag() + x2.imag(); // x[d1] = x1 + x2
		x[d1 + Np2].real() = x1.real() - x2.real();
		x[d1 + Np2].imag() = x1.imag() - x2.imag(); // x[d1 + Np2] = x1 - x2

		return;
	}

	ifft(X, x, Np2, offset, d1, step << 1);
	ifft(X, x, Np2, offset + step, d1 + Np2, step << 1);

	for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

		a = q + d1;
		b = a + Np2;

		x1 = x[a];
		x2.real() = x[b].real() * WN[c].real() + x[b].imag() * WN[c].imag();
		x2.imag() = x[b].imag() * WN[c].real() - x[b].real() * WN[c].imag(); // x2 = x[b]*WN*[c]

		x[a].real() = x1.real() + x2.real();
		x[a].imag() = x1.imag() + x2.imag(); // x[a] = x1 + x2
		x[b].real() = x1.real() - x2.real();
		x[b].imag() = x1.imag() - x2.imag(); // x[b] = x1 - x2
	}

	if (step != 1)
		return;

	_1pN = 1.0 / N;
	for (q = 0; q < N; q++) {
		x[q].real() = x[q].real() * _1pN;
		x[q].imag() = x[q].imag() * _1pN; // x[q] = x[q]/N
	}
}
