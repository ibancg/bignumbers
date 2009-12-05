
#ifndef __BIGNUMBER_H_
#define __BIGNUMBER_H_

// "Arbitrarily" big numbers with fixed-point arithmetic.
//
// A BCD system with sign and modulus representation and decimal adjustment
// was chosen, so each figure represents a decimal digit.

#include "config.h"

class BigNumber {

 public:

  TBC   C[NCIF];
  bool  isPositive; // positive/!negative flag
  
  // Constructors.

  // Creates an empty bignumber
  BigNumber();

  // Creates a bignumber from a string, example: N = BigNumber("-1786.059e36");
  BigNumber(char *);
  
  // Visualization.
  void show();

  // Conversions between bignumbers and floats
  friend void flt2Bn   (FLT flt, BigNumber& bn);
  friend void bn2Flt   (BigNumber& bn, FLT& flt);

  // Tests if two BNs are equal
  friend bool equals (BigNumber& A, BigNumber& B);

  // Compares two BNs and returns the number of coincident digits
  friend int  compare(BigNumber& A, BigNumber& B);
  
  // First non-zero digit index.
  friend int  findFirstNonZeroDigitIndex (BigNumber& A);

  // Copies a BN into another
  friend void copy  (BigNumber& A, BigNumber& B);

  // Operations.

  // Computes C = A + B. If the result is zero, the sign can be explicitly set.
  // A, B and C are overlappables.
  friend void add	(BigNumber& A, BigNumber& B, BigNumber& C, bool sign = true);

  // Computes C = A - B. If the result is zero, the sign can be explicitly set.
  // A, B and C are overlappables.
  friend void sub	(BigNumber& A, BigNumber& B, BigNumber& C, bool sign = true);

  // Computes C = A*B
  friend void mul   (BigNumber& A, BigNumber& B, BigNumber& C);

  // Computes B = A^2
  friend void sqr   (BigNumber& A, BigNumber& B);

  // Computes B = 1/A
  friend void inv   (BigNumber& A, BigNumber& B);

  // Computes C = A/B
  friend void div   (BigNumber& A, BigNumber& B, BigNumber& C);

  // Performs a shift left operation of n digits. The result will be B = A*10^n
  // A and B are overlappables
  friend void shl   (BigNumber& A, BigNumber& B, int n);

  // Performs a shift left operation of n digits. The result will be B = A/10^n
  // A and B are overlappables
  friend void shr   (BigNumber& A, BigNumber& B, int n);

  // Computes the square root
  friend void sqrt  (BigNumber& A, BigNumber& B);

  // Computes the quartic root
  friend void sqrt4 (BigNumber& A, BigNumber& B);
};

#endif
