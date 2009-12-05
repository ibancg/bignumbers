#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <math.h>

#ifndef flt_t
#define flt_t double
#endif

// Complex number.

class Complex {

public:
	flt_t r, i;
};

#endif
