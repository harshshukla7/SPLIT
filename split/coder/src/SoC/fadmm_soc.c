
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
static double resi[nDual];
static double resi_tmp[nDual];

//static double rp[nDual];            // previous step Dual error
//static double sp[nPrimal];          // previous step primal error
static double  Eta, Etainv, iEta, ibeta_k, ibeta_1, beta_k, beta, beta_1, c_nes_prev, c_nes;

#ifdef precond
static double workDual_scale[nDual]; // Tempoary workdual scaled with E^-1
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
    //zero_vector(rp, nDual);
    //zero_vector(sp, nPrimal);
    
}

// Function declaration
void solve(Sol *sol, double par[nParam], const Opt *opt)
{
    
    double rDual, rPrimal, rDual_prev, rPrimal_prev;
    int itr, i;
    
    Eta = 0.999;
    Etainv = 1/Eta;
    beta_k = 1;
    iEta = 1/Eta;
    
    rDual = 10 ;
    rPrimal = 10 ;
    
    
    
    
    loadData();
    // Compute: l = pL*par + l_const, etc
    custom_compute_parametric(l, f, b, par);
    
    // Set kktRHS[nPrimal+1:end] = pB*par + b
    copy_vector(kktRHS+nPrimal, b, nEqCon);
    
    //printf("KKT RHS is %f, %f, %f, and\n b is %f, %f, %f and\n nprimal is %f", kktRHS[0], kktRHS[2], kktRHS[1], b[0], b[1], b[2], nPrimal );
    // Compute termination tolerances
    const double DualTol = (opt->dualTol);
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
    
    
    #ifdef warm_start_x
    copy_vector(x, warm_x, nPrimal);
    #endif
    
    #ifdef warm_start_y
    copy_vector(y, warm_y, nDual);
    copy_vector(y_hat, warm_y, nDual);
    #endif
    
   
    
    #ifdef warm_start_lambda
    copy_vector(lambda, warm_dual_lambda, nDual);
    
    copy_vector(lambda_hat, warm_dual_lambda, nDual);
    #endif
    
    
    
        // c = norm(L*x+l-yd)^2;
    
#ifdef precond
    custom_mult_Ld(workDual, x);
    forall(nDual) resi[i] = workDual[i]+ld[i]-y[i];
#endif
#ifndef precond
    custom_mult_L(workDual, x);
    forall(nDual) resi[i] = workDual[i]+l[i]-y[i];
#endif
    
    
    
    c_nes_prev = sum_squared(resi, nDual);
    c_nes = c_nes_prev ;
    

    
    
    long double start_kkt, end_kkt, sum_time, start_total, end_total;
    sum_time = 0;
    //opt->MAXITR
    start_total = split_tic();
    
    for(itr=0; itr< opt->MAXITR ; itr++)
    {
        
        
        
        rDual_prev = rDual;
        rPrimal_prev = rPrimal;
        
        
        //////////////////////////
        /////////////////////////
        /////////////////////////
        ////////////////////////
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
        forall(nDual) workDual[i] = (rho*(-ld[i] + y_hat[i])) - lambda_hat[i];
        
        // kktRHS[1:nPrimal] = L'*workDual
        custom_mult_Ldtrans(kktRHS, workDual);
#endif
        
#ifndef precond
        forall(nDual) workDual[i] = (rho*(-l[i] + y_hat[i])) - lambda_hat[i];
        
        char *WD_tp = "workdual";
        print_vector(WD_tp, workDual, nDual);
        printf("rho is %f \n ", rho);
        
        custom_mult_Ltrans(kktRHS, workDual);
#endif
        forall(nPrimal) kktRHS[i] -= f[i];
        
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
        // Compute workDual = rho*(-l + y - lambda)
        forall(nDual) workDual[i] = (rho*(-ld[i] + y_hat[i])) - lambda_hat[i];
        
        /////////////////////////////////////
        ///// delete the following lines
        //char *WD_tp = "li";
        //print_vector(WD_tp, y_hat, nDual);
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
        forall(nDual) workDual[i] = ((rho*(-l[i] + y_hat[i])) - lambda_hat[i]);
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
        
        start_kkt = split_tic(); // for timings
        // Solve the KKT system
        custom_solve_kkt(x, kktRHS);
        end_kkt = split_toc(start_kkt);
        //printf("end_KKt is %Lf\n",end_kkt);
        sum_time = end_kkt + sum_time;
        
        
        
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
        forall(nDual) resi[i] = workDual[i] + ld[i];
        
        //forall(nDual) workDual[i] = workDual[i] + ld[i] + (lambda_hat[i]*rhoinv);
        
        forall(nDual) workDual[i] = resi[i] + (lambda_hat[i]*rhoinv);
        forall(nDual) workDual_scale[i] = Einv_vec[i]*workDual[i];
        
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *workdual="Einv_vec";
        //print_vector(before_prox, Einv_vec, (nDual));
        // printf("\n \n \n  first element of workdual is %f Einv is %f and workdual_scale is %f \n \n \n ",workDual[0],Einv_vec[0],workDual_scale[0]);
        
        //////////////// delete part ends here
        
        // Evaluate prox functions y = prox(workDual)
        copy_vector(prev_y, y, nDual);
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
        forall(nDual) resi[i] = workDual[i] + l[i];
        
        
        // workDual = lambda + workDual + l
        
        //forall(nDual) workDual[i] = workDual[i] + l[i] + (lambda_hat[i]*rhoinv);
        
        forall(nDual) workDual[i] = resi[i] + (lambda_hat[i]*rhoinv);
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *workdual="Einv_vec";
        //print_vector(before_prox, Einv_vec, (nDual));
        // printf("\n \n \n  first element of workdual is %f  and \n \n \n ",workDual[0]);
        
        //////////////// delete part ends here
        // Evaluate prox functions y = prox(workDual)
        copy_vector(prev_y, y, nDual);
        custom_prox(y, workDual);
        
        //char *y_after_print="After";
        //print_vector(y_after_print, y, nDual);
        
        //char *y_before_print="Before";
        //print_vector(y_before_print, prev_y, nDual);
        
        
        
#endif
        forall(nDual) resi[i] = resi[i] - y[i];
        
        //char *y_after_print="Resi";
        //print_vector(y_after_print, resi, nDual);
        
        
        forall(nDual) resi_tmp[i] = prev_y[i] - y[i];
        
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
        
        
        copy_vector(prev_lambda, lambda, nDual);
        forall(nDual) lambda[i] = rho*(workDual[i] - y[i]);
        
        
        ////////////////////////////////////
        ////////// delete the following print line
        ////////////////////////////////////////
        //char *lambda_after_print="lambda_after_solve";
        //print_vector(lambda_after_print, lambda, nDual);
        //////////////// delete part ends here
        
        
        /**********************************************************************
         *
         *  // Nestrov's Relaxation [Step 4]
         *
         **********************************************************************/
        
        
        // c = 1/rho * norm(lam-hat.lam)^2 + rho*norm(yd-hat.yd)^2;
        c_nes_prev = c_nes;
        
        forall(nDual) r[i] = y[i] - y_hat[i];
        c_nes = sum_squared(r,nDual);
        c_nes *=rho;
        
        forall(nDual) r[i] = lambda[i] - lambda_hat[i];
        rDual = sum_squared(r,nDual);
        rDual *= rhoinv;
        
        c_nes = c_nes+rDual;
        
        
        //printf("c is %f prev c is %f and eta is %f \n",c_nes,c_nes_prev,Eta);
        
        if (c_nes < Eta*c_nes_prev)
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
            copy_vector(lambda_hat, prev_lambda, nDual);
            copy_vector(y_hat, prev_y, nDual);
            c_nes = Etainv * c_nes_prev;
        }
        
        /**********************************************************************
         *
         *  // Adaptive step-size [Step 5]
         *
         **********************************************************************/
        
#ifdef adaptive
        
        rPrimal = sum_absolute(resi,nDual);
        
#ifdef precond
        custom_mult_residual(s, resi_tmp);
        
#endif
        
#ifndef precond
        custom_mult_Ltrans(s, resi_tmp);
#endif
        forall(nPrimal) s[i] *= rho;
        rDual = sum_absolute(s,nPrimal);
        
        
        
        if (rPrimal > 10*rDual)
        { rho = 2*rho;
          rho_inv = 1/rho;
        }
        
        if (rDual > 10*rPrimal)
        { rho = 0.5*rho;
          rho_inv = 1/rho;
        }
#endif
        
        //char *lambda_after_print="lambda_after_solve";
        //print_vector(lambda_after_print, lambda, nDual);
        // Check convergence
        /**********************************************************************
         *
         *  // Convergence Check [Step 6]
         *
         **********************************************************************/
        if (itr % opt->ITR_PER_CONV_TEST == 0)
        {
#ifndef adaptive
            rPrimal = sum_absolute(resi,nDual);
#endif
            
            if (rPrimal < PrimalTol)
            {
                
                //char *lambda_after_print="workDual_convergence";
                //print_vector(lambda_after_print, resi_tmp, nDual);
                
                
                // rPrimal = ||L'*workDual - rho||^2
#ifndef adaptive
#ifdef precond
                custom_mult_residual(s, resi_tmp);
                
#endif
                
#ifndef precond
                custom_mult_Ltrans(s, resi_tmp);
#endif
                forall(nPrimal) s[i] *= rho;
                rDual = sum_absolute(s,nPrimal);
                //printf("Dual residual is %f",rDual);
#endif
                
                if (rDual < DualTol){
                    
                    printf("The fADMM loop is broken due to PrimalTolSquares\n");
                    break;
                }
            }
        }
    }
    end_total = split_toc(start_total);
    /////////////////////////////////////
    // Copy to solution structure
    copy_vector(sol->primal, x, nPrimal);
    copy_vector(sol->dual, lambda, nDual);
    copy_vector(sol->aux_prim, y, nDual);
    sol->itr = itr;
    sol->rDual = rDual;
    sol->rPrimal = rPrimal;
    sol->time_sol = sum_time;
    sol->time_total = end_total;
    
#ifdef suitesparse_linsolve
    
    free(Li_ss);
    free(Lx_ss);
#endif
}
