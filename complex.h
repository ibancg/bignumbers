
#ifndef _COMPLEX_H_
#define _COMPLEX_H_


#include <math.h>

#ifndef FLT
#define FLT double
#endif

// Aritmetica compleja sencilla.
class CPX {

 public:

   FLT r, i;


   friend void SumaCPX(CPX A, CPX B, CPX &R)
   {
      R.r = A.r + B.r;
      R.i = A.i + B.i;
   }

   friend void RestaCPX(CPX A, CPX B, CPX &R)
   {
      R.r = A.r - B.r;
      R.i = A.i - B.i;
   }

   friend void MulCPX(CPX A, CPX B, CPX &R)
   {
      R.r = A.r*B.r - A.i*B.i;
      R.i = A.i*B.r + A.r*B.i;
   }
};


#endif
