/*
 * Solve parametric convex optimization problem using ADMM
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx 
 *
 * f = pF*par + f_const
 * l = pL*par + l_const
 * b = pB*par + b_const
 *
 * Sizes:
 *  nParam  = length(par)
 *  nPrimal = size(A,2)
 *  nDual   = size(L,1)
 *  nEqCon  = size(A,1)
 *
 */

#ifndef __ADMM__H__
#define __ADMM__H__

#include <stdio.h>
#include <string.h>
#include "user_matrix_ops.h"

#ifndef __probData_h__
#include "user_probData.h"
#endif

#include "splitTimer.h" // needed only in case of timings

// Structure to return solution
struct Sol { 
  double primal[nPrimal]; // Primal variable
  double dual[nDual];     // Dual variable
  int    itr;             // Number of iterations
  double rDual;           // Dual residual
  double rPrimal;         // Primal residual
  long double time_sol;            // to pass measured time
  long double time_total;  // total time taken
  double aux_prim[nDual];  // y
};
typedef struct Sol Sol;

// Options
struct Opt { 
  double primalTol;           // Primal tolerance
  double dualTol;             // Dual tolerance
  unsigned int   MAXITR;              // Max allowed number of iterations
  unsigned int   ITR_PER_CONV_TEST;   // Number of iterations between convergence tests
}; 
typedef struct Opt Opt;
// const Opt Opt_Default = {1e-4, 1e-4, 1000, 10};

// Initialize values of all variables to zero on first call
void initialize();

// Solve the problem
void solve(Sol *sol, double par[nParam], const Opt *opt);

#endif
