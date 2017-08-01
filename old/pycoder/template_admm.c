/*
 * Solve parametric convex optimization problem using ADMM
 *
 [[[cog
 import genSplit as split
 split.gen_header()
 ]]]
 //[[[end]]]
 */

#include "admm.h"

// Variable Definitions 
static double x[nPrimal+nEqCon]; // Extra rows are working space for solving KKT system
static double y[nDual];
static double lambda[nDual];
static double prev_lambda[nDual];
static double prev_y[nDual];

static double workDual[nDual];      // Working memory of size dual
static double kktRHS[nPrimal+nEqCon]; // RHS when solving the KKT system
static double r[nDual];             // Dual error
static double s[nPrimal];           // Primal error


// Load required data into static memory
[[[cog

# If pB is non-zero, then
#    b = pB*par + const_b
# else
#    b = const_b
par_vars = (('l','pL'),('f','pF'),('b','pB'))
for var,parVar in par_vars:
  if split.dat[parVar].size == 0:
    split.load_data(var)
  else:
    split.load_data(var, defName='const_%s' % var)
    split.load_data(parVar)
    numvars = split.dat[var].shape[0]
    cog.outl("static double %s[%i];" % (var,numvars))


# Load all specified constants into static memory
split.load_constants()

# Compute the inverse of the KKT matrix and load it
# TODO : Add some options here where the user can decide which method to use to solve the KKT system

import numpy as np
from numpy.linalg import inv

Q   = split.dat['Q']
rho = split.constants['rho'] # Note: We are assuming rho is constant here
L   = split.dat['L']
A   = split.dat['A']

# K = [Q+rho*L'*L A'; A zeros(size(A,1))];
K = np.vstack((np.hstack((Q+rho*np.dot(L.T,L), A.T)), np.hstack((A, np.zeros((A.shape[0], A.shape[0]))))))
invK = inv(K)
# Set near-zero elements to zeros
invK[np.abs(invK) < 1e-6] = 0
split.load_data('invK', invK)
]]]
//[[[end]]]

// Initialize values of all variables
void zero_vector(double *vec, int len) {for(int i=0; i<len; i++) vec[i]=0.0;}
void initialize()
{
  zero_vector(x, nPrimal);
  zero_vector(y, nDual);
  zero_vector(lambda, nDual);
  zero_vector(prev_lambda, nDual);
  zero_vector(prev_y, nDual);
  zero_vector(workDual, nDual);
}

// Function declaration
[[[cog
split.gen_function_prototype()
]]]
//[[[end]]]
{
  double rDual, rPrimal;
  int itr, i;

  [[[cog
  par_vars = (('l','pL'),('f','pF'),('b','pB'))
  for var,parVar in par_vars:
    if split.dat[parVar].size > 0:
      # This is a parametric variable - compute it
      split.matvec_product(parVar, 'par', var)
      numvars = split.dat[var].shape[0]
      cog.outl("forall(%i) %s[i] += const_%s[i];" % (numvars, var, var))
  ]]]
  [[[end]]]

  // Set kktRHS[nPrimal+1:end] = pB*par + b
  copy_vector(kktRHS+nPrimal, b, nEqCon);

  // print_vector("b", b, nEqCon);
  // print_vector("kktRHS", kktRHS, nPrimal+nEqCon);

	for(itr=0; itr<MAXITR; itr++) 
	{
    copy_vector(prev_lambda, lambda, nDual);
    copy_vector(prev_y, y, nDual);

    // Compute workDual = rho*(-l + y - lambda)
    forall(nDual) workDual[i] = rho*(-l[i] + y[i] - lambda[i]);

    // kktRHS[1:nPrimal] = L'*workDual
    [[[cog
    # Generate code to multiply a sparse constant L by the vector workDual
    split.matvec_product('L', 'workDual', 'kktRHS', transpose=True)
    ]]]
    [[[end]]]

    // print_vector("workDual", workDual, nDual);
    // print_vector("kktRHS", kktRHS, nPrimal+nEqCon);


    forall(nPrimal) kktRHS[i] -= f[i];

    // Solve the KKT system
    // Solve using dense mat-vec product, since we pre-inverted the KKT matrix
    matvec_product(invK, kktRHS, x, nPrimal+nEqCon, nPrimal+nEqCon);

    // print_vector("x", x, nPrimal+nEqCon);

    // workDual = L*x
    [[[cog
    split.matvec_product('L', 'x', 'workDual')
    ]]]
    [[[end]]]

    // workDual = lambda + workDual + l
    forall(nDual) workDual[i] = lambda[i] + workDual[i] + l[i];

    // Evaluate prox functions y = prox(workDual)
    [[[cog
    split.gen_prox('workDual', 'y', 'rho')
    ]]]
    [[[end]]]

    // Dual update
    forall(nDual) lambda[i] = workDual[i] - y[i];

    // Check convergence
    if (itr % ITR_PER_CONV_TEST == 0)
    {
      // Dual tolerance rDual = (lambda-prev_lambda)^2
      forall(nDual) r[i] = lambda[i] - prev_lambda[i];
      rDual = sum_squared(r,nDual);

      if (rDual < DUALTOL*DUALTOL) 
      {
        forall(nDual) workDual[i] = y[i] - prev_y[i];
        [[[cog
        # Generate code to multiply a sparse constant L by the vector workDual
        split.matvec_product('L', 'workDual', 's', transpose=True)
        ]]]
        [[[end]]]
        forall(nPrimal) s[i] *= -rho;
        rPrimal = sum_squared(s,nPrimal);

        if (rPrimal < PRIMALTOL*PRIMALTOL)
          break;
      }
    }
	}

  // Copy to solution structure
  copy_vector(sol->primal, x, nPrimal);
  copy_vector(sol->dual, lambda, nDual);
  sol->itr = itr;
  sol->rDual = rDual;
  sol->rPrimal = rPrimal;
}