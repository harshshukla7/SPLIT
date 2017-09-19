
/*
 * Solve parametric convex optimization problem using ADMM
 *
 * min 0.5 x'*Q*x + f'*x + sum w_i prox_i(y_i)
 * s.t. Ax = b
 *      y  = Lx 
 *
 */

#include "PDA.h"

// Variable Definitions 
static double x[nPrimal]; // Extra rows are working space for solving KKT system
//static double y[nDual];
static double lambda[nDual]; 
static double lambda_tilde[nDual]; 
static double xy_work[nPrimal]; // kkt-RHS, also used as workspace
static double prev_x[nPrimal]; 
static double p[nPrimal];
static double y[nPrimal];
static double x_bar[nPrimal];
static double prev_lambda[nDual]; 
static double workDual[nDual];      // Working memory of size dual
static double dx[nPrimal];
static double dlambda[nDual];
static double b_work[nEqCon];
static double b_temp[nEqCon];

static double r[nPrimal];             // Dual error
static double s[nDual];         // Primal error


//static double prim_eqcons = nPrimal+nEqCon;
//static double dual_eqcons = nDual+nEqCon;

// Initialize values of all variables
void zero_vector(double *vec, int len) {for(int i=0; i<len; i++) vec[i]=0.0;}
void initialize()
{  
  
  
  zero_vector(x, nPrimal);
  zero_vector(xy_work, nPrimal);
  zero_vector(prev_x, nPrimal);
  zero_vector(p, nPrimal);
  zero_vector(y, nPrimal);
  zero_vector(x_bar, nPrimal);
  zero_vector(dx, nPrimal);
  
  zero_vector(lambda, nDual); // May bot be used delete after words
  zero_vector(prev_lambda, nDual); // May bot be used delete after words
  zero_vector(workDual, nDual);
  zero_vector(dlambda, nDual); // May bot be used delete after words
  zero_vector(lambda_tilde, nDual); // May bot be used delete after words
  zero_vector(b_work,nEqCon);
  zero_vector(b_temp,nEqCon);
  
}

// Function declaration
void solve(Sol *sol, double par[nParam], const Opt *opt)
{
    
  double rDual, rPrimal;
  int itr, i;
    
  
  loadData();
  

// Compute: l = pL*par + l_const, etc
  custom_compute_parametric(l, f, b, par);

  // Set xy_work[nPrimal+1:end] = pB*par + b
//  copy_vector(xy_work+nPrimal, b, nEqCon);
  
  //printf("KKT RHS is %f, %f, %f, and\n b is %f, %f, %f and\n nprimal is %f", xy_work[0], xy_work[2], xy_work[1], b[0], b[1], b[2], nPrimal );
  // Compute termination tolerances
  const double DualTolSquared   = (opt->dualTol)*(opt->dualTol);
  const double PrimalTolSquared = (opt->primalTol)*(opt->primalTol);
  const double DualTol = opt->dualTol;
  const double PrimalTol = opt->primalTol;
  
  
   #ifdef warm_start_x
    copy_vector(x, warm_x, nPrimal);
    #endif
    
   // #ifdef warm_start_y
   // copy_vector(y, warm_y, nDual);
   // #endif
    
    #ifdef warm_start_lam
    copy_vector(lambda, warm_lambda, nDual);
    #endif
   
    long double sum_time, start_total, end_total;
    sum_time = 0.0;
    //opt->MAXITR
    start_total = split_tic(); 
    
    
    ////////////////////////////////////////
    //////////////// Warm start
    ////////////////////////////////////////
    
    copy_vector(x, sol->primal, nPrimal+nEqCon);
    copy_vector(lambda, sol->dual, nDual);
    
  //opt->MAXITR
  for(itr=0; itr < opt->MAXITR; itr++) 
  {
   
   // Solve the KKT system
    
  /*************************** Solve KKT system ***********************
    ******************************************************************
   *********************Step 1(a) If adaptive
   *********************Step (1,1) adaptive and precond
   **********************Step (1,2) adaptive and not precond
   *********************Step 1(b) If not adaptive
   *********************Step (1,3) not adaptive and precond
   **********************Step (1,4) not adaptive and not precond  */
    
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step (1,1) adaptive and precond
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    

    // Data required : L, prec.T, dat.f, TQ, TQtau
    // (0) prev_x = x;
    // (i) L'*lambda = xy_work -- mat-vec prod
    // (ii) L'*lambda + dat.f' = xy_work -- vec-vec add
    // (iii) prec.T*(L'*lambda + dat.f') = p -- vec-vec prod since prec.T is diagonal
    // (iv) tau*(TQ*x_prev) = x mat-vec
    // (v) (eye())*x - tau*(TQ*x_prev)  = x -- vec-vec product -substraction  
    // (vi) (TQ- (eye()/tau)*x) + (prec.T*(L'*lambda + dat.f')) = y = y + p -- vec-vec add
    // (vii) tau*(TQ- (eye()/tau)*x) + (prec.T*(L'*lambda + dat.f')) = y = y*tau -- scaler-vec prod
    
    
    #ifdef adaptive
    #ifdef precond
    
    copy_vector(prev_x, x, nPrimal);
    custom_mult_Ltrans(xy_work, lambda);
    forall(nPrimal) xy_work[i] += f[i];
    forall(nPrimal) y[i] = prec_T[i]*xy_work[i];
    custom_mult_TQ(x, prev_x);
    forall(nPrimal) x[i] = prev_x[i] - (x[i]*tau);
    forall(nPrimal) y[i] = x[i] - (y[i]*tau);
    
    #endif
    #endif
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step (1,2) adaptive and not precond
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    // Data required : L, dat.f, dat.Q, Qtau
    // (0) prev_x = x;
    // (i) L'*lambda = y -- mat-vec prod
    // (ii) L'*lambda + dat.f' = y -- vec-vec add
    // (iii) (Q- (eye()/tau) = Qtau do this smartly since its diagonal --vec substra 
    // (v) (Q- (eye()/tau))*x = p = Qtau*x-- mat-vec product  
    // (vi) (- Q + (eye()/tau)*x) + ((L'*lambda + dat.f')) = y = y + p  -- vec-vec add
    // (vii) tau*(Q- (eye()/tau)*x) + ((L'*lambda + dat.f')) -- scaler-vec prod
    
    
    #ifdef adaptive
    #ifndef precond
    
    copy_vector(prev_x, x, nPrimal);
    custom_mult_Ltrans(xy_work, lambda);
    forall(nPrimal) xy_work[i] += f[i];
    custom_mult_Q(x, prev_x);
    forall(nPrimal) x[i] =  prev_x[i] - (x[i]*tau);
    forall(nPrimal) y[i] = x[i] - (xy_work*tau);
    
    #endif
    #endif
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step (1,3) not adaptive and precond
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    // Data required : L, (tauT), dat.f, ITQ
    // (0) prev_x = x;
    // (i) L'*lambda = xy_work -- mat-vec prod
    // (ii) L'*lambda + dat.f' = xy_work -- vec-vec add
    // (iii) TauT*(L'*lambda + dat.f') = p = TauT*y -- mat-vec prod
    // (iv) (ITQ)*x = y -- mat-vec product  
    // (v) (ITQ*x) + (tauT*(L'*lambda + dat.f')) = y = y + p -- vec-vec add
    
    #ifndef adaptive
    #ifdef precond
    
    copy_vector(prev_x, x, nPrimal);
    custom_mult_Ltrans(xy_work, lambda);
    forall(nPrimal) xy_work[i] += f[i];
    forall(nPrimal) y[i] = tau_T[i]*xy_work[i];
    custom_mult_ITQ(x, prev_x);
    forall(nPrimal) y[i] = x[i] - (y[i]);
    
    #endif
    #endif
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step (1,4)  not adaptive and not precond
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    // Data required : L, dat.f, ItauQ
    // (0) prev_x = x;
    // (i) L'*lambda = y -- mat-vec prod
    // (ii) L'*lambda + dat.f' = y -- vec-vec add
    // (iii) (ItauQ)*x = p -- mat-vec product  
    // (iv) (ItauQ)*x) + ((L'*lambda + dat.f')) = y = y + p -- vec-vec add
    // (vii) tau*(Q- (eye()/tau)*x) + ((L'*lambda + dat.f')) y = rho*y-- scaler-vec prod
    
    #ifndef adaptive
    #ifndef precond
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *x_at_start="x_at_Start";
    //print_vector(x_at_start, x, (nPrimal));
    //////////////// delete part ends here
    
    copy_vector(prev_x, x, nPrimal);
    custom_mult_Ltrans(xy_work, lambda);
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *x_aft_copy="x_aft_copy";
    //print_vector(x_aft_copy, x, (nPrimal));
    //////////////// delete part ends here
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *kkt_RHS="KKT_RHS";
    //print_vector(kkt_RHS, xy_work, (nPrimal));
    //////////////// delete part ends here
    
    
    forall(nPrimal) xy_work[i] += f[i];
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *kkt_RHS2="KKT_RHS2";
    //print_vector(kkt_RHS2, xy_work, (nPrimal));
    //////////////// delete part ends here
    
    
    custom_mult_ItauQ(x, prev_x);
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *x_before="x_bef";
    //print_vector(x_before, x, (nPrimal));
    //////////////// delete part ends here
    //////////////// delete part ends here
    
    forall(nPrimal) y[i] = x[i] - (xy_work[i]*tau);
    
    //printf(" \n \n \n tau is %f \n \n \n",tau);
    #endif
    #endif
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step 2 Dynamics update
    ////////////////// #ifdef dynamics_update 
    ////////////////// step(2,1) datA*y = p
    ////////////////// step(2,2) datb-datA*y = p = p - dat.b
    ////////////////// step(2,3)  p = y + invA'AA'*(dat.b-dat.A*y) = y + invA'AA'*p
    ////////////////// #endif
    ////////////////// #ifndef
    ////////////////// p = y; 
    ////////////////// #endif
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *y_before_dynamics="y_before_dynamics";
    //print_vector(y_before_dynamics, y, (nPrimal));
    //
    //////////////// delete part ends here
        
    
    #ifdef dynamics_update
    
    custom_mult_A(b_work, y);
    forall(nEqCon) b_work[i] = b[i] - b_work[i];
    
    custom_chol_solve(b_temp, b_work);
    custom_mult_Atrans(p,b_temp);
    forall(nPrimal) p[i] = y[i] + p[i];
     //custom_mult_AAA(xy_work, b_work);     
    //        forall(nPrimal) p[i] = y[i] + xy_work[i];
    #endif
    
    #ifndef dynamics_update
    copy_vector( p, y,nPrimal);
    #endif
    
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *p_after_dynamics="p_after_dynamics";
    //print_vector(p_after_dynamics, p, (nPrimal));
    //////////////// delete part ends here

    
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step 3 Update
    ////// (i) x_bar = 2*p-x_prev -- vector substraction
    ////// (ii) p = (p-x) -- vec substraction
    /////  (ii) x = x + theta*p  -- scaler-vector addition and vector addition
    
    forall(nPrimal) x_bar[i] = 2*p[i] - prev_x[i] ;
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *xbar_after_dynamics="xbar_bar";
    //print_vector(xbar_after_dynamics, x_bar, (nPrimal));
    //////////////// delete part ends here
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *x_before_final="xbefore_final";
    //print_vector(x_before_final, x, (nPrimal));
    //////////////// delete part ends here

    
    //printf("\n\n\n\n thetha is %f \n\n\n\n",theta);
    forall(nPrimal) x[i] = prev_x[i] + theta*(p[i]-prev_x[i]);
        
    
    ////////////////////////////////////
    ////////// delete the following print line 
    ////////////////////////////////////////
    //char *x_after_print="x_after_solve";
    //print_vector(x_after_print, x, (nPrimal));
    //////////////// delete part ends here

    
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step 4 Prox. Step
    ////// (i) lambda_prev = lambda;
    /////  (ii) p = L*x_bar mat-vec
    ////   (iii) p =  p + l
    ///// ifdef precond
    ////    (iv) lambda_tilde = prec.P*p
    ////    (v) lambda_tilde  = lambda + rho*lambda_tilde
    ////  endif
    ///    ifndef precond
    /////  (iv) lambda_tilde = lambda +  rho*p - vector addi _ scaler addition 
    ////   (v)  compute prox: mud = conjprox (...)
    ////   (vi) lambda = lambda + theta(mud(ind)'-lambd)
    
    copy_vector(prev_lambda, lambda, nDual);
    custom_mult_L(workDual, x_bar);
    
    #ifdef precond
    forall(nDual) lambda_tilde[i] = lambda[i] + rho*(prec_P[i]*(workDual[i]+l[i])); 
    
    #endif
    
    #ifndef precond
    forall(nDual) lambda_tilde[i] = lambda[i] + rho*(workDual[i]+l[i]);
    #endif
    
    
    //////////////////////////////////////////////////////
    ///////////////////////// Delete the following lines 
    
    
    //char *lambda_tilde_hs="lambda_tilde";
    //print_vector(lambda_tilde_hs, lambda_tilde, nDual);
    
    
    ///////////////////////////
    
    //printf("\n \n \n \n So the values of 46 to 50 are %f %f %f %f %f \n \n \n", lambda_tilde[44], lambda_tilde[45],lambda_tilde[46],lambda_tilde[47],lambda_tilde[48],lambda_tilde[49]);
    custom_prox(workDual, lambda_tilde);
    
    //printf("\n \n \n \n So the values of workdual  46 to 50 are %f %f %f %f %f \n \n \n", workDual[44], workDual[45],workDual[46],workDual[47],workDual[48],workDual[49]);
    forall(nDual) lambda[i] = lambda[i] + theta*(workDual[i]-lambda[i]);
    
    
    //////////////////////////////////////////////////////
    ///////////////////////// Delete the following lines 
    //char *WD_prox="wd_prox";
    //print_vector(WD_prox, workDual, nDual);
    
    //char *lambda_after_print="lambda_after_solve";
    //print_vector(lambda_after_print, lambda, nDual);
    
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step 5 Convergence check
    ///////// (i)  step 1: dx = x - prev.x -- vector substraction
    ////////  (ii) step 2: dlambd = lambda - prev_lambda -- vector substraction
    
    //////    (iii) step 3: r = (DinvL)*dlamd
    //////    (iv)  step 4: s = (QtauD)*dx
    //////     (v)  step 5: r = r + s
    /////      (vi) step 6: rPrimal = norm(r,1)
    
    /////////  (vii) step 7: s = (EL)*dx
    ////////   #ifdef adaptive
    ///////    (viii)step 8(a): y = (Einv)*dlamd
    //////          step 8(b): y = rho*dlambd 
    //////     #endif
    /////      #ifndef adaptive
    //////     (viii)step 8: y = (rhoEinv)*dlambd
    //////     #endif
    //////     (ix) step 9: s = (EL)*dx
    /////      (x) step 10: rDual = norm(s)
    
    forall(nPrimal) dx[i] = x[i] - prev_x[i];
    forall(nDual) dlambda[i] = lambda[i] - prev_lambda[i];
    
    //////////////////////////////////////////////////////
    ///////////////////////// Delete the following lines 
    
    
    //char *dxx="dx";
    //print_vector(dxx, dx, nPrimal);
    
    //char *dlam="dlam";
    //print_vector(dlam, dlambda, nDual);
    
    
    // finding r if adaptive
    
    #ifdef adaptive
    
    custom_mult_Q(r, dx);
    
    #ifdef precond
    
    forall(nDual) r[i] = r[i] - (prec_D[i]*(tau_inv)*dx[i]);
            
    #endif
    
    #ifndef precond
    
    forall(nDual) r[i] = r[i] - (dx[i]*tau_inv);
    
    #endif
    
    #endif 
    
    // finding r if not adaptive
    
    #ifndef adaptive
    
    custom_mult_QTD(r,dx);
    
    #endif
    
    custom_mult_Ltrans(y, dlambda);
    
    
    #ifdef precond
    
    forall(nPrimal) y[i] = prec_Dinv[i]*y[i];
    
    #endif
    
    forall(nPrimal) r[i] = r[i] + y[i];
    
    
    // finding s
    
     custom_mult_L(workDual, dx);
    
     #ifdef precond
     
     forall(nDual) workDual[i] = prec_E[i]*workDual[i];
     forall(nDual) s[i] = -rho_inv*(dlambda[i]/prec_E[i]);
     
     #endif
     
     
     #ifndef precond
     
     forall(nDual) s[i] = -rho_inv*(dlambda[i]);
     
     #endif
     
    // char *sip="si_print";
    //print_vector(sip, s, nDual);
    
    //char *wdp="WD";
    //print_vector(wdp, workDual, nDual);
     
     
    forall(nDual) s[i] = s[i] + workDual[i];
     
    
     rPrimal = sum_absolute(r,nPrimal);
     rDual = sum_absolute(s,nDual);
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////// Step 6 Adaptive step size
     
     #ifdef adaptive
     
     #ifdef Lipschitz
     prev_tau = tau;
     prev_alpha = alpha;
     prev_rho = rho;
     
     tau = tau /(1-alpha);
     rho = rho *(1-alpha);
     
     alpha = alpha * eta;
     
     tau = tau * (1-alpha);
     rho = rho /(1-alpha);
     
     
     alpha *= eta;
     
     if(tau*tau*rho*rho >= ((svd_QT*tau+eps)/sqrt_Mu)*((svd_QT*tau+eps)/sqrt_Mu)){
         
         tau = prev_tau;
         rho = prev_rho;
         alpha = prev_alpha
                 
     }
     
     #endif
     
     if((rPrimal < (s_adapt*rDual)/Delta) ){
      
         tau = tau * (1-alpha);
         rho = rho / (1-alpha) ;
         alpha = alpha * eta;
     }
     
     
     if((rPrimal >= (s_adapt*rDual)/Delta) & (rPrimal <= (s_adapt*rDual)*Delta) ){
      
         tau = tau;
         rho = rho;
         alpha = alpha;
     }
         
     
     #endif
  
    // Check convergence
    if (itr % opt->ITR_PER_CONV_TEST == 0)
    {
      // Dual tolerance rDual = (lambda-prev_lambda)^2
      
      if (rDual < DualTol) 
      {
        
        if (rPrimal < PrimalTol){
            
            printf("The loop is broken due to PrimalTolSquares\n");
            break;
        }
      }
    }
  }
  end_total = split_toc(start_total);
  
  // Copy to solution structure
  copy_vector(sol->primal, x, nPrimal);
  copy_vector(sol->dual, lambda, nDual);
  copy_vector(sol->aux_prim, y, nPrimal);
  sol->itr = itr;
  sol->rDual = rDual;
  sol->rPrimal = rPrimal; 
  sol->time_total = end_total;
  sol->time_sol = 0.0; // Only when measuring linear solve time ! 
}
