
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// conversiones   FLT <-> BigNumber.

#define NUM_DIGITOS 20
#define MAX_EXP     200


void FLT2BN(FLT x, BigNumber &X)
{
  static char buffer[40];
  char        *s = buffer;
  
  sprintf(s, "%20.20f", x);
  memset(X.C, 0, NCIF*sizeof(TBC));
  
  X.positivo = (*s != '-');
  
  unsigned int ls;
  int          i, pp;

  if (!X.positivo) s++;
    
  char *e = strchr(s, 'e');
  ls = ((e) ? (e - s) : strlen(s));

  char *bp;
  int exponente = ((e) ? strtol(&e[1], &bp, 0) : 0);

  if ((exponente > MAX_EXP) || (exponente < -MAX_EXP)) {
    printf("ERROR (FLT2BN): exponente fuera de rango.\n");
    exit(255);
  }

  char *punto = strchr(s, '.');
  pp = ((punto) ? (punto - s) : ls); // posición del punto. 
  
  for (i = 0; i < pp; i++)
    X.C[NFRC + pp - i - 1 + exponente] = s[i] - 48;
  
  if (pp != ls) for (i = pp + 1; i < ls; i++)
    X.C[NFRC - (i - pp) + exponente] = s[i] - 48;
}

//-----------------------------------------------

void BN2FLT(BigNumber &X, FLT &x)
{
  register int i;
  int          n;

  /* Si el orden de magnitud del número es n, empiezo a iterar en 10^(n/2)*/
  for (i = NCIF - 1, n = -1; (i >= 0) && (n == -1); i--)
    if (X.C[i]) n = i;
  
  if (n == -1) {
    x = 0.0;
    return;
  }

  x = 0.0;
  i = n - NUM_DIGITOS + 1;
  if (i < 0) i = 0;

  double exponente = pow(10, i - NFRC);

  for (; i <= n; i++) {
    x += X.C[i]*exponente;
    exponente *= 10.0;
  }
}
