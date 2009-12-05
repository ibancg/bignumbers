
#include "config.h"
#include "complex.h"
#include "FFT.h"


/*
  TRANSFORMADAS RÁPIDAS DE FOURIER.
  ---------------------------------

  (ñ) 2000 Ibán Cereijo Graña.

  Última revisión : 30-03-2001

  Nota : No se hace uso de funciones de aritmética compleja. Basta con que
  la clase de complejos (aquí llamada CPX), tenga dos variables PUBLICAS
  llamadas r e i (que corresponden lógicamente a la parte real y a la parte
  imaginaria)
 */

// factores de fase.
CPX WN[NCIF];

/*
  Optimización:
  -------------

  Crea la tabla de los factores de fase.
  WN[i] = exp(-j*k)       con k = 0..pi  (N muestras)
*/

void CreaWN()
{
  FLT alpha;

  for (int i = 0; i < NCIF; i++) {
    alpha = -i*M_PI/NCIF;
    WN[i].r = cos(alpha);
    WN[i].i = sin(alpha);
  }
}

/* -------------------------------
   transformada rápida de Fourier.
   ------------------------------- */
void FFT(CPX *x, CPX *X, ULI N, ULI offset, ULI d1, ULI step)
{
  CPX           X1, X2;
  ULI           Np2 = (N >> 1); // N/2
  register ULI  a, b, c, q;

  if (N == 2) { // Mariposa para N = 2;
    
    X1 = x[offset];
    X2 = x[offset + step];
    X[d1].r       = X1.r + X2.r;
    X[d1].i       = X1.i + X2.i; // X[d1] = X1 + X2
    X[d1 + Np2].r = X1.r - X2.r;
    X[d1 + Np2].i = X1.i - X2.i; // X[d1 + Np2] = X1 - X2
    return;
  }

  FFT(x, X, Np2, offset, d1, step << 1);
  FFT(x, X, Np2, offset + step, d1 + Np2, step << 1);
    
  for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

    a = q + d1;
    b = a + Np2;

    X1 = X[a];
    X2.r = X[b].r*WN[c].r - X[b].i*WN[c].i;
    X2.i = X[b].r*WN[c].i + X[b].i*WN[c].r; // X2 = X[b]*WN[c]

    X[a].r = X1.r + X2.r;
    X[a].i = X1.i + X2.i; // X[a] = X1 + X2
    X[b].r = X1.r - X2.r;
    X[b].i = X1.i - X2.i; // X[b] = X1 - X2
  }
}

/* ---------------------------------------
   transformada rápida inversa de Fourier.
   --------------------------------------- */
void IFFT(CPX *X, CPX *x, ULI N, ULI offset, ULI d1, ULI step)
{
  CPX           x1, x2;
  ULI           Np2 = (N >> 1); // N/2
  FLT           _1pN = 1.0/N;
  register ULI  a, b, c, q;

  if (N == 2) { // Mariposa para N = 2;
    
    x1 = X[offset];
    x2 = X[offset + step];
    x[d1].r       = x1.r + x2.r;
    x[d1].i       = x1.i + x2.i; // x[d1] = x1 + x2
    x[d1 + Np2].r = x1.r - x2.r;
    x[d1 + Np2].i = x1.i - x2.i; // x[d1 + Np2] = x1 - x2

    return;
  }

  IFFT(X, x, Np2, offset, d1, step << 1);
  IFFT(X, x, Np2, offset + step, d1 + Np2, step << 1);
    
  for (q = 0, c = 0; q < (N >> 1); q++, c += step) {

    a = q + d1;
    b = a + Np2;

    x1 = x[a];
    x2.r = x[b].r*WN[c].r + x[b].i*WN[c].i;
    x2.i = x[b].i*WN[c].r - x[b].r*WN[c].i; // x2 = x[b]*WN*[c] (OJO al conjugado!!!)

    x[a].r = x1.r + x2.r;
    x[a].i = x1.i + x2.i; // x[a] = x1 + x2
    x[b].r = x1.r - x2.r;
    x[b].i = x1.i - x2.i; // x[b] = x1 - x2
  }

  if (step != 1) return;

  _1pN = 1.0/N;
  for (q = 0; q < N; q++) {
    x[q].r = x[q].r*_1pN;
    x[q].i = x[q].i*_1pN; // x[q] = x[q]/N
  }
}
