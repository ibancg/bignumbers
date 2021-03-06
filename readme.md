
DESCRIPTION
-----------

Library designed for deal with arbitrary-precision arithmetic. Also, a simple
example is provided.

DETAILS
-------

The operations implemented in this version are:

 - Addition and subtraction.
 - Multiplication
 - Division
 - Square root.
 - Quartic root.
 
The number of significant digits must be set beforehand in the properties
BigNumber::N_DIGITS and BigNumber::N_FRAC_DIGITS.

In the file config.h you will find the following switches for testing the
performance of the different algorithms. When all of them are switched on, the
performance is maximum:

 - FFT_MUL_ALGORITHM : FFT multiplication algorithm.
 - INVERSE_DIV_ALGORITHM : Division through multiplication by inverse.
 - INVERSE_NEWTON_SQRT_ALGORITHM : Inverse square algorithm
 - INVERSE_NEWTON_SQRT4_ALGORITHM : inverse quartic root algorithm.


COMPILING
---------

Compile it entering the command "make". You will need the GCC compiler and the
GNU Make tool.


RUNNING THE TEST PROGRAM
------------------------

	Usage: bignumber

As an application example, the program evaluates over one million decimals of PI
using a fourth order Borwein's algorithm (1985). The result, with 1038090
decimals, will be dumped to the file 'output.txt'. The last 7 decimals are
incorrect due to precision issues (1038083 digits are correct).

For carrying out the computation, 2^20 N_DIGITS are needed. It can take you 
near one hour to perform the calculation in a normal computer (2GHz CPU). You
can try with less digits. For example, you can compute near 32768 (2^15)
decimals of PI in about 15 seconds.
