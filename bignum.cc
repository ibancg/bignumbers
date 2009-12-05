
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "bignum.h"

// Constructores.
BigNumber::BigNumber()
{
}

// ej: N = BigNumber("-1786.059e36");
BigNumber::BigNumber(char *s)
{
  memset(C, 0, NCIF*sizeof(TBC));

  positivo = (*s != '-');

  unsigned int ls;
  int          i, pp;

  if (!positivo) s++;

  char *e = strchr(s, 'e');
  ls = ((e) ? (e - s) : strlen(s));

  char *bp;
  int exponente = ((e) ? strtol(&e[1], &bp, 0) : 0);

  char *punto = strchr(s, '.');
  pp = ((punto) ? (punto - s) : ls); // posición del punto. 
  
  for (i = 0; i < pp; i++)
    C[NFRC + pp - i - 1 + exponente] = s[i] - 48;
  
  if (pp != ls) for (i = pp + 1; i < ls; i++)
    C[NFRC - (i - pp) + exponente] = s[i] - 48;
}

//------------------------------- VISUALIZACION ------------------------------

/*
   Suponemos que cada cifra interna del formato codifica una cifra decimal,
 (poca optimización de memoria).
*/
void BigNumber::Mostrar()
{
  int                i, j;
  unsigned short int nc = 0;
  unsigned long int  c, ch;
  bool               z = true;

  if (!positivo) printf("-");
   
  for (i = NCIF - 1; i >= NFRC; i--) {
    
    if (z && (C[i])) z = false;
    if (!z) {
      printf("%c", C[i] + 48);
      nc++;
    }
  }
  
  if (z) { // caso especial: el número es el 0.
    printf("0");
    nc++;
  }
  
  for (j = 0; j < NFRC; j++) if (C[j]) break;

  // parte decimal.
  if (j != NFRC) printf(".");
  for (i = NFRC - 1; i >= j; i--) {
    printf("%c", C[i] + 48);
    nc++;
  }
  
  printf("::(%u cifras)\n", nc);
}

//------------------------------------------------------------------------

// SHL en potencias de 10.
// A y B solapables.
void ShlBN(BigNumber &A, BigNumber &B, int d)
{
  if (d >= 0) {
    memcpy(B.C + d, A.C, NCIF - d*sizeof(TBC));
    memset(B.C, 0, d*sizeof(TBC));
  } else {
    memcpy(B.C, A.C - d, NCIF + d*sizeof(TBC));
    memset(&B.C[NCIF + d], 0, -d*sizeof(TBC));
  }
  
  B.positivo = A.positivo;	  
}

// SHR en potencias de 10.
// A y B solapables.
void ShrBN(BigNumber &A, BigNumber &B, int d)
{
  ShlBN(A, B, -d);	
}

// Transferir.
void TraBN(BigNumber &A, BigNumber &B)
{
  memcpy(B.C, A.C, NCIF*sizeof(TBC));
  B.positivo = A.positivo;
}

//------------------------------------------------------------------------

// devuelve el índice de la primera cifra no nula (empezando por arriba).
int IndicePC(BigNumber &A)
{
  for (register int i = NCIF - 1; i >= 0; i--)
    if (A.C[i]) return i;

  return -1; // es un 0
}

bool ComparaBN(BigNumber &A, BigNumber &B)
{
  for (register int i = 0; i < NCIF; i++)
    if (A.C[i] != B.C[i]) return false;  
  return true;
}

// devuelve el número de cifras que coinciden.
int Compara2BN(BigNumber &A, BigNumber &B)
{
  for (register int i = NCIF - 1; i >= 0; i--)
    if (A.C[i] != B.C[i]) return (NCIF - i);

  return NCIF;
}
