
#ifndef __BIGNUMBER_H_
#define __BIGNUMBER_H_

/*
--------------------------------------------------------------------------------
   Formato num�rico y aritm�tica asociada de n�meros "arbitrariamente" grandes
 con convenio de punto fijo. 
 
   Se ha escogido un sistema BCD con ajuste decimal en cada d�gito interno,
 as� cada cifra interna codificar� 1 d�gito decimal. La aritm�tica binaria
 natural se rechaza por las complicaciones que presenta a la hora de visua-
 lizar el n�mero (algoritmo de orden n^2).
 
   En cuanto al signo se usa el sistema signo y m�dulo.

--------------------------------------------------------------------------------   
 (�) 2000 Ib�n Cereijo Gra�a
--------------------------------------------------------------------------------   
*/

#include "config.h"

class BigNumber {

 public:

  TBC   C[NCIF];
  bool  positivo; // flag positivo/!negativo
  
  // Constructores.
  BigNumber();
  BigNumber(char *);
  
  // Visualizaci�n.
  void Mostrar();

  // conversiones
  friend void FLT2BN   (FLT, BigNumber &);
  friend void BN2FLT   (BigNumber &, FLT &);

  friend bool ComparaBN (BigNumber &, BigNumber &);
  friend int  Compara2BN(BigNumber &, BigNumber &);
  
  // �ndice de la primera cifra no nula.
  friend int  IndicePC (BigNumber &);

  // Operaciones.
  friend void SumaBN   (BigNumber &, BigNumber &, BigNumber &, bool = true);
  friend void RestaBN  (BigNumber &, BigNumber &, BigNumber &, bool = true);
  friend void MulBN    (BigNumber &, BigNumber &, BigNumber &);
  friend void SqrBN    (BigNumber &, BigNumber &);
  friend void InvBN    (BigNumber &, BigNumber &);
  friend void DivBN    (BigNumber &, BigNumber &, BigNumber &);
  friend void ShlBN    (BigNumber &, BigNumber &, int);
  friend void ShrBN    (BigNumber &, BigNumber &, int);
  friend void TraBN    (BigNumber &, BigNumber &);
  friend void SqrtBN   (BigNumber &, BigNumber &);
  friend void Sqrt4BN  (BigNumber &, BigNumber &);
};

#endif
