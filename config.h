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

#ifndef __CONFIG_H_
#define __CONFIG_H_

// Configuration file.

#include <math.h>

// debug mode.
#define DEBUG

typedef unsigned char bcd_t; // BCD digit type.

// advanced configuration (algorithm switches).

#define FFT_MUL_ALGORITHM
#define INVERSE_DIV_ALGORITHM
#define INVERSE_NEWTON_SQRT_ALGORITHM
#define INVERSE_NEWTON_SQRT4_ALGORITHM

#endif
