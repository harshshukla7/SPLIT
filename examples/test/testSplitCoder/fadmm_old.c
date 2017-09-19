
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
static double lambda_hat[nDual];
static double prev_y[nDual];
static double y_hat[nDual];

static double workDual[nDual];      // Working memory of size dual
static double kktRHS[nPrimal+nEqCon]; // RHS when solving the KKT system
static double r[nDual];             // Dual error
static double s[nPrimal];           // Primal error

//static double rp[nDual];            // previous step Dual error
//static double sp[nPrimal];          // previous step primal error
static double beta, beta_k, Ep, ibeta_k, beta_1, ibeta_1;

// Initialize values of all variables
void zero_vector(double *vec, int len) {for(int i=0; i<len; i++) vec[i]=0.0;}
void initialize()
{  
  
  
  zero_vector(x, nPrimal);
  zero_vector(y, nDual);
  zero_vector(lambda, nDual);
  zero_vector(prev_lambda, nDual);
  zero_vector(lambda_hat, nDual);
  zero_vector(prev_y, nDual);
  zero_vector(y_hat, nDual);
  zero_vector(workDual, nDual);
  //zero_vector(rp, nDual);
  //zero_vector(sp, nPrimal);
  
}

// Function declaration
void solve(Sol *sol, double par[nParam], const Opt *opt)
{
    
  double rDual, rPrimal, rDual_prev, rPrimal_prev;
  beta = 0;
  Ep = 1;
  beta_k = 1;
  rDual = 0 ;
  rPrimal = 0 ; 

 int itr, i;
    
  
  loadData();
  // Compute: l = pL*par + l_const, etc
  custom_compute_parametric(l, f, b, par);

  // Set kktRHS[nPrimal+1:end] = pB*par + b
  copy_vector(kktRHS+nPrimal, b, nEqCon);
  
  //printf("KKT RHS is %f, %f, %f, and\n b is %f, %f, %f and\n nprimal is %f", kktRHS[0], kktRHS[2], kktRHS[1], b[0], b[1], b[2], nPrimal );
  // Compute termination tolerances
  const double DualTolSquared   = (opt->dualTol)*(opt->dualTol);
  const double PrimalTolSquared = (opt->primalTol)*(opt->primalTol);
  
  
  //for (int i=0; i<28; i++){
  //  printf("Ap_ss matrix entry %d is %f \n",i ,Ax_ss );
//}
  
  custom_compute_prefactor();
  /*for (int i=0; i<28; i++){
    printf("D_ss before in the main  just out side prefactor %d is %f\n", i, D_ss[i]);
}
   *
   *check here
   *
   *
for(int  i= 0; i < N_ss; i++){
printf("Lnz[%d] is %d \n", i, Lnz[i]);
}

for (int i=0; i<38; i++){
    printf("Lx_ss before in the main  just out side prefactor %d is %f\n", i, Lx_ss[i]);
} */
  //for (int i=0; i<28; i++){
  //  printf("Diagonal matrix entry %d is %f \n",i ,D_ss );
//}
  
  
  //opt->MAXITR
  for(itr=0; itr< opt->MAXITR; itr++) 
  {
   
    
    
    rDual_prev = rDual;
    rPrimal_prev = rPrimal;

    // Compute workDual = rho*(-l + y - lambda)
    forall(nDual) workDual[i] = rho*(-l[i] + y_hat[i] - lambda_hat[i]);

    // kktRHS[1:nPrimal] = L'*workDual
    custom_mult_Ltrans(kktRHS, workDual);

    forall(nPrimal) kktRHS[i] -= f[i];
    
   
   
    // Solve the KKT system
    



    custom_solve_kkt(x, kktRHS);
    
   
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    char *x_after_print="x_after_solve";
    print_vector(x_after_print, x, (nPrimal));
    //////////////// delete part ends here

    // workDual = L*x
    custom_mult_L(workDual, x);

    


    // workDual = lambda + workDual + l
    forall(nDual) workDual[i] = lambda_hat[i] + workDual[i] + l[i];

    // Evaluate prox functions y = prox(workDual)
    copy_vector(prev_y, y, nDual);
    
    custom_prox(y, workDual);
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    char *y_after_print="y_after_solve";
    print_vector(y_after_print, y, nDual);
    //////////////// delete part ends here

    // Dual update
    
    copy_vector(prev_lambda, lambda, nDual);
    forall(nDual) lambda[i] = +workDual[i] - y[i];
     
    
    // Nestrov's Relaxation
    // Dual tolerance rDual = (lambda-prev_lambda)^2
    forall(nDual) r[i] = lambda[i] - prev_lambda[i];
    rDual = sum_squared(r,nDual);
    
    forall(nDual) workDual[i] = y[i] - prev_y[i];

        // rPrimal = ||L'*workDual - rho||^2
        custom_mult_Ltrans(s, workDual);
        forall(nPrimal) s[i] *= -rho;
        rPrimal = sum_squared(s,nPrimal);
      
    Ep = fmax(rPrimal_prev, rDual_prev) - fmax(rPrimal, rDual);
    
    if (Ep > 0)
    {
        beta = beta_k;
        beta_k = (1 + sqrt(4*beta*beta+1))*0.5;
        ibeta_k = 1/beta_k;
        beta_1 = (beta-1);
        ibeta_1 = ibeta_k * beta_1;
        
        forall(nDual) workDual[i] = lambda_hat[i] + workDual[i] + l[i];
        forall(nDual) y_hat[i] = y[i] + (ibeta_1*(y[i]-prev_y[i]));
        forall(nDual) lambda_hat[i] = lambda[i] + (ibeta_1 * (lambda[i]-prev_lambda[i]));
        //y_hat = y + (beta_1)*(y - prev_y)*ibeta_k;
        //lambda_hat = lambda + (beta_1)*(lambda - prev_lambda)*ibeta_k;
    }
    else 
    {
        beta_k = 1;
        copy_vector(lambda_hat, lambda, nDual);
        copy_vector(y_hat, y, nDual);
        
    }
        
        
    //char *lambda_after_print="lambda_after_solve";
    //print_vector(lambda_after_print, lambda, nDual);
    // Check convergence
    if (itr % opt->ITR_PER_CONV_TEST == 0)
    {
      
      

      if (rDual < DualTolSquared) 
      {
        

        if (rPrimal < PrimalTolSquared){
            
            printf("The fADMM loop is broken due to PrimalTolSquares\n");
            break;
        }
      }
    }
  }

  // Copy to solution structure
  copy_vector(sol->primal, x, nPrimal);
  copy_vector(sol->dual, lambda, nDual);
  sol->itr = itr;
  sol->rDual = rDual;
  sol->rPrimal = rPrimal;
  
  // Harsh Edits here.. deallocate the memory
  
  free(Li_ss);
  free(Lx_ss);
  
}
