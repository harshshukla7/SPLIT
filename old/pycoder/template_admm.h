/*
 * Solve parametric convex optimization problem using ADMM
 *
 [[[cog
 import genSplit as split
 split.gen_header()
 ]]]
 //[[[end]]]
 */

#ifndef __ADMM__H__
#define __ADMM__H__

#include <stdio.h>
#include <string.h>
#include "matrix_ops.h"

// Variable sizes
[[[cog
(r,c) = split.dat['pB'].shape
cog.outl('#define nParam  %i' % c)

(r,c) = split.dat['L'].shape
cog.outl('#define nPrimal %i' % c)
cog.outl('#define nDual   %i' % r)

(r,c) = split.dat['A'].shape
cog.outl('#define nEqCon  %i' % r)
]]]
[[[end]]]

// Structure to return solution
struct Sol { 
  double primal[nPrimal]; // Primal variable
  double dual[nDual];     // Dual variable
  int    itr;             // Number of iterations
  double rDual;           // Dual residual
  double rPrimal;         // Primal residual
};
typedef struct Sol Sol;

// Initialize values of all variables
void initialize();

// Function declaration
///[[[cog
split.gen_function_prototype()
cog.outl(';')
///]]]
//[[[end]]]


#endif
