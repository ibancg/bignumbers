
#include <string.h>
#include <stdio.h>
#include <iostream>

#include "bignum.h"

using namespace std;

// Cálculo del inverso por el método de Newton.
void InvBN(BigNumber &A, BigNumber &B)
{
  static BigNumber x1, x2;
  static BigNumber DOS("2");
  int              i;
  FLT              d;
  int              ipc = 0; // índice primera cifra.
  bool             stop;

# ifdef DEBUG
  cout << "INV";
  cout.flush();
#endif

  B.positivo = A.positivo;
  memset(B.C, 0, NCIF*sizeof(TBC)); // limpiamos C.

  // cuento el número de cifras enteras.
  for (i = NCIF - 1; (i >= 0) && !ipc; i--) if (A.C[i]) ipc = i;
  B.C[NFRC - (ipc - NFRC) - 1] = 1;
  
  for (int k = 0;; k++) {
    
    MulBN(A, B, x1);
    RestaBN(DOS, x1, x2);
    MulBN(B, x2, x1);

# ifdef DEBUG
    cout << '.';
    cout.flush();
#endif

    if (ComparaBN(B, x1)) break;

    TraBN(x1, B);
  }

# ifdef DEBUG
  cout << endl;
#endif
}


#ifdef ALGORITMO_DIV_INVERSO

void DivBN(BigNumber &A, BigNumber &B, BigNumber &C)
{
  InvBN(B, C);
  MulBN(A, C, C);
}

#else
//-------------------------------------------------------------

// El método del cole. solapable A con B o C.
void DivBN(BigNumber &A, BigNumber &B, BigNumber &C)
{
  int                i, j, n;
  int                NCA = 0, NCB = 0, NCAB;
  static TBC         AA[NCIF + NFRC];


  // AA va a ser modificado.
  memcpy(&AA[NFRC], A.C, NCIF*sizeof(TBC)); // copio A desplazado en AA.
  memset(AA, 0, NFRC*sizeof(TBC));

  for (i = NCIF - 1; i >= 0; i--) {
    if ((!NCA) && (A.C[i])) NCA = i;
    if ((!NCB) && (B.C[i])) NCB = i;
  }

  NCA += NFRC; // AA est  desplazado NFRC cifras respecto a A.

  if (!NCB && !B.C[0]) {
    printf("ERROR: división por 0\n");
    exit(255);
  }
  
  memset(C.C, 0, NCIF*sizeof(TBC));
  C.positivo = !(A.positivo ^ B.positivo);

  if (NCB > NCA) return;
  
  NCAB = NCA - NCB;
  
  bool               menor;
  char               r, c;
 
  i = 1;  
  for (i = 0; i <= NCAB; i++) {
    
    // Compruebo si AA es menor que B desplazado.
    for (n = 0;; n++) {
      
      menor = false;
      for (j = NCB + 1; j >= 0; j--) {
	
	if (AA[j + NCAB - i] == B.C[j]) continue;
	if (AA[j + NCAB - i] <  B.C[j]) menor = true;
	break;
      }
      if (menor) break;
      
      // le resto B desplazado a AA.
      for (c = 0, j = 0; j <= NCB + 1; j++) {
		  
        r = AA[j + NCAB - i] - (B.C[j] + c);
        c = (r < 0) ? 1 : 0; 
        AA[j + NCAB - i] = (r + 10*c);  // r % 10
      }
    }
	
    C.C[NCAB - i] = n;
  }		  
}
#endif
