
/*
 * Solve parametric convex optimization problem using ADMM
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx
 *
 */


#include "user_ama.h"

//#include "splitTimer.h"
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
static double r[nDual];             // Primal error
static double s[nPrimal];           // Dual error

static double beta, beta_k, Ep, ibeta_k, beta_1, ibeta_1;

#ifdef precond
static double workDual_scale[nDual]; // Tempoary workdual scaled with E^-1
#endif

#ifdef adaptive_restart
double ad_rest ;
#endif

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
    
#ifdef precond
    zero_vector(workDual_scale, nDual);
#endif
    
    
}

// Function declaration
void solve(Sol *sol, double par[nParam], const Opt *opt)
{
    
    double rDual, rPrimal, rDual_prev, rPrimal_prev;
    int itr, i;
    
    Ep = 1;
    beta_k = 1; //current_beta
    beta = 1; // previous_beta
    rDual = 10 ;
    rPrimal = 10 ;
    
 //   loadData();
    // Compute: l = pL*par + l_const, etc
#ifndef precond
    custom_compute_parametric(l, f, b, par);
#endif
    
#ifdef precond
    custom_compute_parametric(ld, f, b, par);
#endif
    
    // Set kktRHS[nPrimal+1:end] = pB*par + b
    copy_vector(kktRHS+nPrimal, b, nEqCon);
    
    //printf("KKT RHS is %f, %f, %f, and\n b is %f, %f, %f and\n nprimal is %f", kktRHS[0], kktRHS[2], kktRHS[1], b[0], b[1], b[2], nPrimal );
    // Compute termination tolerances
    //const double DualTolSquared   = (opt->dualTol)*(opt->dualTol);
    //const double PrimalTolSquared = (opt->primalTol)*(opt->primalTol);
    
    const double DualTol   = (opt->dualTol);
    const double PrimalTol = (opt->primalTol);
    
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
    
    long double start_kkt, end_kkt, sum_time, start_total, end_total; 
    sum_time = 0;    
    //opt->MAXITR
  //  start_total = split_tic();
    for(itr=0; itr < opt->MAXITR ; itr++)
    {
        // Step 1: Nesterov relaxation
//
        
        beta_k = (1 + sqrt(4*beta*beta+1))*0.5;
        beta = beta_k;
        ibeta_k = 1/beta_k;
        beta_1 = (beta-1);
        ibeta_1 = ibeta_k * beta_1;
        
        
        
        forall(nDual) lambda_hat[i] = lambda[i] + (ibeta_1 * (lambda[i]-prev_lambda[i]));
        
        //  char *lambda_hat_init_after_print="first_lambda_hat";
        //     print_vector(lambda_hat_init_after_print, lambda_hat, (nDual));
        
        
        /**********************************************************************
         *
         *  Solve KKT for adaptive with precondition case [Step 2]
         *
         **********************************************************************/
        
        // Compute workDual = rho*(-l + y - lambda)
        forall(nDual) workDual[i] = (-lambda_hat[i]);
        
#ifdef precond
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ldtrans(kktRHS, workDual);
#endif
        
#ifndef precond
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ltrans(kktRHS, workDual);
#endif
        
        
        forall(nPrimal) kktRHS[i] -= f[i];
        // Set kktRHS[nPrimal+1:end] = pB*par + b
        copy_vector(kktRHS+nPrimal, b, nEqCon);
    
        
        
        // Solve the KKT system
        
        
        
        //printf("KKT Solve");
        //start_kkt = split_tic();
        custom_solve_kkt(x, kktRHS);
        //end_kkt = split_toc(start_kkt); 
        //printf("end_KKt is %Lf\n",end_kkt);
        //sum_time = end_kkt + sum_time;
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
         //char *x_after_print="x_after_solve";
         //print_vector(x_after_print, x, (nPrimal));
        //////////////// delete part ends here
        
        /**********************************************************************
         *
         *  Prox Computation [Step 3]
         *
         **********************************************************************/
        
        /**********************************************************************
         *
         * Compute Prox for with precondition case [Substep (2,1)]
         *
         ***********************************************************************/
        
        copy_vector(prev_y, y, nDual);
        
#ifdef precond
// workDual = L*x
        custom_mult_Ld(workDual, x);
        
        forall(nDual) r[i] = workDual[i] + ld[i];
        // workDual = lambda + workDual + l
        forall(nDual) workDual[i] = (rhoinv*lambda_hat[i]) + r[i];
        
        forall(nDual) workDual_scale[i] = Einv_vec[i]*workDual[i];
        
        // Evaluate prox functions y = prox(workDual)
        custom_prox(y, workDual_scale);
        
        forall(nDual) y[i] = E_vec[i]*y[i];
        
        forall(nDual) r[i] = r[i] - y[i]; // just use the new y for computing residual
#endif
        
        
#ifndef precond
        // workDual = L*x
        custom_mult_L(workDual, x);
        
        forall(nDual) r[i] = workDual[i] + l[i]; // first save it here to use later for residual
        // workDual = lambda + workDual + l
        forall(nDual) workDual[i] = (lambda_hat[i]*rhoinv) + r[i];
        
        //char *workDual_print="Workdual";
        //print_vector(workDual_print, workDual, nDual);
        
        // Evaluate prox functions y = prox(workDual)
        custom_prox(y, workDual);
        
         forall(nDual) r[i] = r[i] - y[i]; // just use the new y for computing residual
        
#endif
        
        //printf("rho is %f \n",rho);
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *y_after_print="y_after_solve";
        //print_vector(y_after_print, y, nDual);
        //////////////// delete part ends here
        
        /**********************************************************************
         *
         *  Lagrange multiplier update [Step 4]
         *
         **********************************************************************/
        
        copy_vector(prev_lambda, lambda, nDual);
        // Dual update
        
        // printf("rho is %f \n", rho);
        //char *workDual_print="Workdual";
        //print_vector(workDual_print, workDual, nDual);
        
        
        forall(nDual) lambda[i] = rho*(workDual[i] - y[i]);
        
        //char *lambda_after_print="lambda_after_solve";
        //print_vector(lambda_after_print, lambda, nDual);
        
        /**********************************************************************
         *
         * Adaptive restart [Step 5]
         *
         **********************************************************************/
        //# if 0
#ifdef adaptive_restart
        ad_rest = 0;
        forall(nDual) ad_rest += ( lambda_hat[i] - lambda[i])*(lambda[i] - prev_lambda[i]);
        //printf("adrest is %f",ad_rest);
        if (ad_rest > 0){
            
            copy_vector(lambda_hat, prev_lambda, nDual);
            
            //      char *lambda_after_print="x_after_solve";
            //  print_vector(lambda_after_print, lambda_hat, (nDual));
            
            
            // Compute workDual = rho*(-l + y - lambda)
            forall(nDual) workDual[i] = (-lambda_hat[i]);
            
#ifdef precond
            // kktRHS[1:nPrimal] = L'*workDual
            custom_mult_Ldtrans(kktRHS, workDual);
#endif
            
#ifndef precond
            // kktRHS[1:nPrimal] = L'*workDual
            custom_mult_Ltrans(kktRHS, workDual);
#endif
            
            
            forall(nPrimal) kktRHS[i] -= f[i];
            copy_vector(kktRHS+nPrimal, b, nEqCon);
    
            
            
            // Solve the KKT system
        //    start_kkt = split_tic();
        custom_solve_kkt(x, kktRHS);
       // end_kkt = split_toc(start_kkt); 
        //sum_time = end_kkt + sum_time;
        
        // Next is Prox step
            
#ifdef precond
            //workDual = L*x
            custom_mult_Ld(workDual, x);
            
            forall(nDual) r[i] = workDual[i] + ld[i];
            
            // workDual = lambda + workDual + l
            forall(nDual) workDual[i] = (rhoinv*lambda_hat[i]) + r[i];
            
            forall(nDual) workDual_scale[i] = Einv_vec[i]*workDual[i];
            
            // Evaluate prox functions y = prox(workDual)
            custom_prox(y, workDual_scale);
            
            forall(nDual) y[i] = E_vec[i]*y[i];
            
            forall(nDual) r[i] = r[i] - y[i];
#endif
            
            
#ifndef precond
            // workDual = L*x
            custom_mult_L(workDual, x);
            
            forall(nDual) r[i] = workDual[i] + l[i];
            // workDual = lambda + workDual + l
            forall(nDual) workDual[i] = (rhoinv*lambda_hat[i]) + r[i];
            
            // Evaluate prox functions y = prox(workDual)
            custom_prox(y, workDual);
            
            forall(nDual) r[i] = r[i] - y[i];
            
#endif
            
            copy_vector(prev_lambda, lambda_hat, nDual);
            // Dual update
            forall(nDual) lambda[i] = rho*(workDual[i] - y[i]);
        }
        
#endif
        
//#endif
        /**********************************************************************
         *
         *  Convergence check [Step 6]
         *
         **********************************************************************/
        
        // Check convergence
        if (itr % opt->ITR_PER_CONV_TEST == 0)
        {   rDual_prev = rDual;
            rPrimal_prev = rPrimal;
            
            
            // Dual tolerance rDual = (lambda-prev_lambda)^2
            
            
            rDual = sum_absolute(r,nDual);
            //printf("rDual is %f \n",rDual);
            
            if (rDual < DualTol)
            {
                
                forall(nDual) workDual[i] =  prev_lambda[i] - lambda[i] ;
                
                
#ifdef precond
                custom_mult_residual(s, workDual);
#endif
                
#ifndef precond
                custom_mult_Ltrans(s, workDual);
#endif
                
                rPrimal = sum_absolute(s,nPrimal);
                // printf("rPrimal is %f \n",rPrimal);
                
                
                if (rPrimal < PrimalTol){
                    
                    printf("The loop is broken due to Primal and Dual residual is less than tolerance values.\n");
                    break;
                }
            }
        }
    }
    
    //end_total = split_toc(start_total);
    //////////////////////////////////////
    
    // Copy to solution structure
    copy_vector(sol->primal, x, nPrimal);
    copy_vector(sol->dual, lambda, nDual);
    copy_vector(sol->aux_prim, y, nDual);
    sol->itr = itr;
    sol->rDual = rDual;
    sol->rPrimal = rPrimal;
    sol->time_sol = sum_time;
    sol->time_total = end_total;
    // Deallocate the memory in case of suitesparse
#ifdef suitesparse_linsolve
    free(Li_ss);
    free(Lx_ss);
    
#endif
    
    
}
