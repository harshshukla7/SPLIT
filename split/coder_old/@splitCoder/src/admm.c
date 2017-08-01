/*
 * Solve parametric convex optimization problem using ADMM
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx 
 *
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
void solve(Sol *sol, const double par[nParam], const Opt *opt)
{
  double rDual, rPrimal;
  int itr, i;

  // Compute: l = pL*par + l_const, etc
  custom_compute_parametric(l, f, b, par);

  // Set kktRHS[nPrimal+1:end] = pB*par + b
  copy_vector(kktRHS+nPrimal, b, nEqCon);

  // Compute termination tolerances
  const double DualTolSquared   = (opt->dualTol)*(opt->dualTol);
  const double PrimalTolSquared = (opt->primalTol)*(opt->primalTol);

  for(itr=0; itr<1; itr++) 
  {
    copy_vector(prev_lambda, lambda, nDual);
    copy_vector(prev_y, y, nDual);

    // Compute workDual = rho*(-l + y - lambda)
    forall(nDual) workDual[i] = rho*(-l[i] + y[i] - lambda[i]);

    // kktRHS[1:nPrimal] = L'*workDual
    custom_mult_Ltrans(kktRHS, workDual);

    forall(nPrimal) kktRHS[i] -= f[i];

    // Solve the KKT system
    custom_solve_kkt(x, kktRHS);

    // workDual = L*x
    custom_mult_L(workDual, x);

    // workDual = lambda + workDual + l
    forall(nDual) workDual[i] = lambda[i] + workDual[i] + l[i];

    // Evaluate prox functions y = prox(workDual)
    custom_prox(y, workDual);

    // Dual update
    forall(nDual) lambda[i] = workDual[i] - y[i];

    // Check convergence
    if (itr % opt->ITR_PER_CONV_TEST == 0)
    {
      // Dual tolerance rDual = (lambda-prev_lambda)^2
      forall(nDual) r[i] = lambda[i] - prev_lambda[i];
      rDual = sum_squared(r,nDual);

      if (rDual < DualTolSquared) 
      {
        forall(nDual) workDual[i] = y[i] - prev_y[i];

        // rPrimal = ||L'*workDual - rho||^2
        custom_mult_Ltrans(s, workDual);
        forall(nPrimal) s[i] *= -rho;
        rPrimal = sum_squared(s,nPrimal);

        if (rPrimal < PrimalTolSquared)
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