/*
 * Solve parametric convex optimization problem using AMA
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx 
 *
 * f = pF*par + f_const
 * l = pL*par + l_const
 *
 * Solve parametric convex optimization problem using AMA
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

#ifndef __AMA__H__
#define __AMA__H__

#include <stdio.h>
#include <string.h>
#include "user_matrix_ops.h"
#ifndef __probData_h__
#include "user_probData.h"
#endif

#include "user_foo_data.h"
#include "foo_function_wrapped.h"

// Structure to return solution
struct Sol { 
  data_t_primal_out primal[nPrimal]; // Primal variable
  data_t_primal_out dual[nDual];     // Dual variable
  int    itr;             // Number of iterations
  data_t_primal_out rDual;           // Dual residual
  data_t_primal_out rPrimal;         // Primal residual
  data_t_primal_out aux_prim[nDual];  // y
};
typedef struct Sol Sol;

// Options
struct Opt { 
  data_t_primal_out primalTol;           // Primal tolerance
  data_t_primal_out dualTol;             // Dual tolerance
  int MAXITR;              // Max allowed number of iterations
  int ITR_PER_CONV_TEST;   // Number of iterations between convergence tests
}; 
typedef struct Opt Opt;

// const Opt Opt_Default = {1e-4, 1e-4, 1000, 10};

// Initialize values of all variables to zero on first call

void initialize();

// Solve the problem

//void solve(Sol *sol, double par[nParam], const Opt *opt, data_t_tol_iterates_in tol_iterates_in_int[TOL_ITERATES_IN_LENGTH], data_t_iterates_out iterates_out_int[ITERATES_OUT_LENGTH] );
void solve(Sol *sol, data_t_primal_out par[nParam], const Opt *opt );


#endif

