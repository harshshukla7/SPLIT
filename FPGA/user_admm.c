
/*
 * Solve parametric convex optimization problem using ADMM
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx
 *
 */

#include "user_admm.h"
//#include "ama.h"

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

#ifdef precond
static double workDual_scale[nDual]; // Tempoary workdual scaled with E^-1
#endif


// Initialize values of all variables
void zero_vector(double *vec, int len) {for(int i=0; i<len; i++) vec[i]=0.0;}
void initialize()
{
    zero_vector(x, nPrimal+nEqCon);
    zero_vector(y, nDual);
    zero_vector(lambda, nDual);
    zero_vector(prev_lambda, nDual);
    zero_vector(prev_y, nDual);
    zero_vector(workDual, nDual);
}

// Function declaration
void solve(Sol *sol, double par[nParam], const Opt *opt)
{
    
    double rDual, rPrimal;
    int itr, i;
    
    loadData();
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
     * printf("D_ss before in the main  just out side prefactor %d is %f\n", i, D_ss[i]);
     * }
     *
     *check here
     *
     *
     * for(int  i= 0; i < N_ss; i++){
     * printf("Lnz[%d] is %d \n", i, Lnz[i]);
     * }
     *
     * for (int i=0; i<38; i++){
     * printf("Lx_ss before in the main  just out side prefactor %d is %f\n", i, Lx_ss[i]);
     * } */
    //for (int i=0; i<28; i++){
    //  printf("Diagonal matrix entry %d is %f \n",i ,D_ss );
    //}
    
#ifdef adaptive
    double rho = rho_tmp;
    double rhoinv = rhoinv_tmp;
#endif
   
    
    /* Iterative Loop for solving
     *
     * It should have the following structure:
     *
     *  Step 1 : KKT Solve
     *          Substep (1,1): If adaptive
     *                      subsub step (1,1,1): if preconditioned
     *                      subsub step (1,1,2): if not preconditioned
     *          Substep (1,2): If not adaptive
     *                      subsub step (1,2,1): if preconditioned
     *                      subsub step (1,2,2): if not preconditioned
     *
     *  Step 2 : Prox computation
     *          Substep (2,1): If preconditioning (must be diagonal)
     *          Substep (2,2): If not preconitioned
     *
     *  Step 3: Dual Update
     *
     *  Step 4: \rho update for the case of adaptive
     *
     *  Step 5: Convergence check
     *
     *
     *
     *opt->MAXITR
     */
    #ifdef warm_start_x
    copy_vector(x, warm_x, nPrimal);
    #endif
    
    #ifdef warm_start_y
    copy_vector(y, warm_y, nDual);
    #endif
    
    #ifdef warm_start_lambda
    copy_vector(lambda, warm_dual_lambda, nDual);
    #endif
    
    long double start_kkt, end_kkt, sum_time, start_total, end_total;
    sum_time = 0.0;
    //opt->MAXITR
//    start_total = split_tic();
    for(itr=0; itr < opt->MAXITR ; itr++)
    {
        //printf("\n\n\n\n\n\n\n\n iteration %d \n\n\n\n\n\n\n\n" ,itr);
        copy_vector(prev_lambda, lambda, nDual);
        copy_vector(prev_y, y, nDual);
        
        
        /**********************************************************************
         *
         *  Solve KKT [Step 1]
         *
         **********************************************************************/
        
        
#ifdef adaptive
        
        /**********************************************************************
         *
         *  Solve KKT for adaptive with precondition case [Substep (1,1,1)]
         *
         **********************************************************************/
        
        // step 1: update the rho1 or rho2 based on case I or case II
        
        
        
        // step 1: compute the RHS for KKT system based on new step size
        // step 2: compute updated x
        // 	substep 1 : update division vector smartly depending on the case I or case II
        // 	substep 2 : scale the matrix with division
        // 	substep 3 : matrix-vector multiplication
        // 	substep 4 : matrix-vector multiplication
        
        
        
#ifdef precond
        forall(nDual) workDual[i] = (rho*(-ld[i] + y[i])) - lambda[i];
        
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ldtrans(kktRHS, workDual);
#endif
        
#ifndef precond
        forall(nDual) workDual[i] = (rho*(-l[i] + y[i])) - lambda[i];
        
        char *WD_tp = "workdual";
        print_vector(WD_tp, workDual, nDual);
        printf("rho is %f \n ", rho);
        
        custom_mult_Ltrans(kktRHS, workDual);
#endif
        forall(nPrimal) kktRHS[i] -= f[i];
        copy_vector(kktRHS+nPrimal, b, nEqCon);
        
        // Print the KKT RHS
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        char *KKT_RHS_beforesovle="KKT_RHS_beforesovle";
        print_vector(KKT_RHS_beforesovle, kktRHS, (nPrimal+nEqCon));
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        //////////////// delete part ends here
        
        custom_mult_XTrans(x,kktRHS);
        
        
        // Print the KKT RHS
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        char *Xtra_rhs="Xtrans_times_RHS";
        print_vector(Xtra_rhs, x, (nPrimal+nEqCon));
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        //////////////// delete part ends here
        
        
        // Print the KKT RHS
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *Xtran_print="Xtrans_vec";
        //print_vector(Xtran_print, XTrans_adapt_vec_p_ss, 86);
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        //////////////// delete part ends here
        
        
#ifdef adap_case_1
        forall(nPrimal) kktRHS[i] = x[i]/(D1_vec[i] + rho);
#endif
        
        
#ifdef adap_case_2
        forall(nPrimal) kktRHS[i] = x[i]/((rho*D1_vec[i]) + 1);
#endif
        
        // Print the KKT RHS
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *divXtra_rhs="divXtrans_times_RHS";
        //print_vector(divXtra_rhs, x, (nPrimal+nEqCon));
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        //////////////// delete part ends here
        
        
        // Print the KKT RHS
        
        custom_mult_X(x,kktRHS);
        
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *step1_final="step1_final_solution";
        //print_vector(step1_final, x, (nPrimal+nEqCon));
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        //////////////// delete part ends here
        
#endif
        
        
#ifndef adaptive
        
        /**********************************************************************
         *
         *  Solve KKT for non adaptive with precondition case [Substep (1,2,1)]
         *
         **********************************************************************/
        
#ifdef precond
        //printf(" \n \n %f is rho. \n It goes in the not adaptive precondition condition \n",rho);
        // Compute workDual = rho*(-l + y - lambda)
        forall(nDual) workDual[i] = (rho*(-ld[i] + y[i])) - lambda[i];
        
        /////////////////////////////////////
        ///// delete the following lines
        //char *WD_tp = "li";
        //print_vector(WD_tp, workDual, nDual);
        //printf("rho is %f \n ", rho);
        //////////////////////////////////////
        
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ldtrans(kktRHS, workDual);
        //
        forall(nPrimal) kktRHS[i] -= f[i];
        copy_vector(kktRHS+nPrimal, b, nEqCon);
        
        // Solve the KKT system``
        custom_solve_kkt(x, kktRHS);
#endif
        
        /**********************************************************************
         *
         * Solve KKT for non adaptive with no precondition case [Substep (1,2,2)]
         *
         ***********************************************************************/
        
#ifndef precond
        
        // Compute workDual = rho*(-l + y - lambda)
        forall(nDual) workDual[i] = ((rho*(-l[i] + y[i])) - lambda[i]);
        
        //printf("\n\n\n rho is %f \n\n\n",rho);
        //char *WD_tp = "workdual";
        //print_vector(WD_tp, workDual, nDual);
        
        
        
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ltrans(kktRHS, workDual);
        //
        forall(nPrimal) kktRHS[i] -= f[i];
        copy_vector(kktRHS+nPrimal, b, nEqCon);
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *KKT_RHS_beforesovle="KKT_RHS_beforesovle";
        //print_vector(KKT_RHS_beforesovle, kktRHS, (nPrimal+nEqCon));
        //////////////// delete part ends here
        
       // start_kkt = split_tic(); // for timings
        // Solve the KKT system
        custom_solve_kkt(x, kktRHS);
        //end_kkt = split_toc(start_kkt);
        //printf("end_KKt is %Lf\n",end_kkt);
        //sum_time = end_kkt + sum_time;
        
        
#endif
        
#endif
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *x_after_print="x_after_solve";
        //print_vector(x_after_print, x, (nPrimal));
        
        //////////////// delete part ends here
        
        
        
        
        
        
        /**********************************************************************
         *
         *  Prox Computation [Step 2]
         *
         **********************************************************************/
        
        /**********************************************************************
         *
         * Compute Prox for with precondition case [Substep (2,1)]
         *
         ***********************************************************************/
#ifdef precond
        // workDual = L*x
        custom_mult_Ld(workDual, x);
        
        
        
        // workDual = lambda + workDual + l
        forall(nDual) r[i] = workDual[i] + ld[i] - y[i];
        
        
        
        // following alpha corresponds to Lx_relax. See boyd for more information
        
#ifndef alpha1 // if alpha is not one
        forall(nDual) workDual[i] = (r[i])*(alpha) + y[i] + lambda[i]*rhoinv;
        //printf("\n \n \n r[0] is %f, alpha is %f, y[0] is %f, lambda[0] is %f, rhoinv is %f  and workdual beofre is %f \n \n \n ",r[0], alpha, y[0], lambda[0], rhoinv,  workDual[0]);
        
#endif
        
#ifdef alpha1 // if alpha is one then do not mulptiply
        forall(nDual) workDual[i] = r[i] + y[i] + lambda[i]*rhoinv;
#endif
        
        forall(nDual) workDual_scale[i] = Einv_vec[i]*workDual[i];
        
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *workdual="Einv_vec";
        //print_vector(before_prox, Einv_vec, (nDual));
        //printf("\n \n \n  first element of workdual is %f Einv is %f and workdual_scale is %f \n \n \n ",workDual[0],Einv_vec[0],workDual_scale[0]);
        
        //////////////// delete part ends here
        
        // Evaluate prox functions y = prox(workDual)
        custom_prox(y, workDual_scale);
        
        forall(nDual) y[i] = E_vec[i]*y[i];
        
#endif
        
        
        /**********************************************************************
         *
         * Compute Prox for with no precondition case [Substep (2,2)]
         *
         ***********************************************************************/
        
        
#ifndef precond
        // workDual = L*x
        custom_mult_L(workDual, x);
        forall(nDual) r[i] = workDual[i] + l[i] - y[i];
        // workDual = lambda + workDual + l
        
        
#ifndef alpha1
        forall(nDual) workDual[i] = (r[i])*(alpha) + y[i] + (lambda[i]*rhoinv);
        
#endif
        
#ifdef alpha1
        forall(nDual) workDual[i] = r[i] + y[i] + lambda[i]*rhoinv;
#endif
        
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *workdual="Einv_vec";
        //print_vector(before_prox, Einv_vec, (nDual));
        //printf("\n \n \n  first element of workdual is %f  and \n \n \n ",workDual[0]);
        
        //////////////// delete part ends here
        // Evaluate prox functions y = prox(workDual)
        custom_prox(y, workDual);
        
        
        
#endif
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *y_after_print="y_after_solve";
        //print_vector(y_after_print, y, nDual);
        //////////////// delete part ends here
        
        /**********************************************************************
         *
         *  Dual Update [Step 3]
         *
         **********************************************************************/
        
        // Dual update
#ifdef gamma1
        forall(nDual) lambda[i] = rho*(+workDual[i] - y[i]);
#endif
        
#ifndef gamma1
        forall(nDual) lambda[i] = gamma*rho*(+workDual[i] - y[i]);
      
        //printf("\n\n\n gamma is %f and rho is %f \n\n\n",gamma,rho);
#endif
        
        
        //char *lambda_after_print="lambda_after_solve";
        //print_vector(lambda_after_print, lambda, nDual);
        
        // Check convergence
        
        /**********************************************************************
         *
         *  Update \rho for the adaptive case [Step 4]
         *
         **********************************************************************/
        
#ifdef adaptive_everystep
        
        // Primal Residual
        forall(nDual) r[i] = r[i] + prev_y[i] - y[i];
        
        
#ifdef precond
        r[i] = Einv_vec[i]*r[i];
#endif
        
        rPrimal = sum_absolute(r,nDual);
        
        // Dual Residual
        
        forall(nDual) workDual[i] = y[i] - prev_y[i];
        
#ifdef precond
        custom_mult_residual(s, workDual);
#endif
        
#ifndef precond
        custom_mult_Ltrans(s, workDual);
#endif
        forall(nPrimal) s[i] *= -rho;
        rDual = sum_absolute(s,nPrimal);
        
        if(rPrimal > c_adapt_step*rDual )
        {rho = 2*rho;
         rhoinv = 1/rho;
        }
        else if (rDual > c_adapt_step*rPrimal ){
            rho = 0.5*rho;
            rhoinv = 1/rho;
        }
        
        if (rPrimal < PrimalTol && rDual < DualTol){
            
            printf("The loop is broken due to residual conditions satisfaction. \n");
            break;
        }
#endif
        
        if (itr % opt->ITR_PER_CONV_TEST == 0)
        {
#ifdef adaptive
            
#ifndef adaptive_everystep
            

            // Primal  Residual
            
            forall(nDual) r[i] = r[i] + prev_y[i] - y[i];
            
#ifdef precond
            r[i] = Einv_vec[i]*r[i];
#endif
            
            rPrimal = sum_absolute(r,nDual);
            
            // Dual Residual
            forall(nDual) workDual[i] = y[i] - prev_y[i];
            
            // rPrimal = ||L'*workDual - rho||^2
#ifdef precond
            custom_mult_residual(s, workDual);
#endif
            
#ifndef precond
            custom_mult_Ltrans(s, workDual);
#endif
            
            forall(nPrimal) s[i] *= -rho;
            
            
            rDual = sum_absolute(s,nPrimal);
            
            
            ////////////////////////////////////
            ////////// delete the following print line
            ////////////////////////////////////////
            //char *rprimal_vec="Primaal_residual_vector";
            //print_vector(rprimal_vec, s, nPrimal);
            
            //printf("\n\n\n\n Primal residual is %f \n\n\n\n",rPrimal);
            
            //char *rdual_vector="Dual_residual_vector";
            //print_vector(rdual_vector, r, nDual);
            
            //printf("\n\n\n\n Dual residual is %f \n\n\n\n",rDual);
            
            
            //char *y_dual_R="y_after_solve";
            //print_vector(y_dual_R, y, nDual);
            //////////////// delete part ends here
            
            
            if(rPrimal > c_adapt_step*rDual )
            {rho = 2*rho;
             rhoinv = 1/rho;
            }
            else if (rDual > c_adapt_step*rPrimal ){
                rho = 0.5*rho;
                rhoinv = 1/rho;
            }
            else{
            }
            
            if (rPrimal < PrimalTol && rDual < DualTol){
                
                printf("The loop is broken due to residual conditions satisfaction. \n");
                break;
            }
            //Primal Residual
#endif
#endif
            
#ifndef adaptive
            
            
            // Primal Residual 
            forall(nDual) r[i] = r[i] + prev_y[i] - y[i];
            
#ifdef precond
            r[i] = Einv_vec[i]*r[i];
#endif
            ////////////////////////////////////
            ////////// delete the following print line
            ////////////////////////////////////////
            // char *r_dual_R="ri_after_solve";
            //print_vector(r_dual_R, r, nDual);
            //char *prev_y_dual_R="prev_y_after_solve";
            //print_vector(prev_y_dual_R, prev_y, nDual);
            //char *y_dual_R="y_after_solve";
            //print_vector(y_dual_R, y, nDual);
            //////////////// delete part ends here
            
            // delete the  following lines
            
            
            
            //
            rPrimal = sum_absolute(r,nDual);
            
            //printf("Primal residual is %f \n",rPrimal);
            
            
            
            
            
            if (rPrimal < PrimalTol)
            {    
                
                // If primal Residual is less than tolerance then checn Dual
                forall(nDual) workDual[i] =  prev_y[i] - y[i];
                
                // rPrimal = ||L'*workDual - rho||^2
#ifdef precond
                custom_mult_residual(s, workDual);
                
#endif
                
#ifndef precond
                custom_mult_Ltrans(s, workDual);
#endif
                forall(nPrimal) s[i] *= rho;
                rDual = sum_absolute(s,nPrimal);
                //printf("Dual residual is %f",rDual);
                if (rDual < DualTol){
                    
                    printf("The loop is broken due to Primal and Dual residuals are less than tolerance values.\n");
                    break;
                }
                
            }
#endif
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
    // deallocate the memory
    //char *y_after_print="y_after_solve";
    //print_vector(y_after_print, y, nDual);
    //printf("Further sum time and total time is %Lf and %Lf and in solution it is %Lf and %Lf",sum_time,end_total,sol->time_sol,sol->time_total);

#ifdef suitesparse_linsolve
    
    free(Li_ss);
    free(Lx_ss);
    
#endif
    
}



