#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// header for mex
#include "mex.h"

// define macros for solvers
// #define ADMM
// #define AMA 
// #define FAMA
// #define FADMM
// #define CPI
#define CPII

// c files for ADMM + header for data
#ifdef ADMM
    #include "solver_admm.c"  
    #include "solve_KKT_admm.c"
    #include "ADMM_matrices.h"
#endif    
        
// c files for AMA + header for data
#ifdef AMA
    #include "solver_ama.c"
    #include "solve_KKT_ama.c"
    #include "AMA_matrices.h"
#endif   
        
// c files for FAMA + header for data
#ifdef FAMA
    #include "solver_fama.c"
    #include "solve_KKT_fama.c"
    #include "FAMA_matrices.h"
#endif   
        
// c files for FADMM + header for data
#ifdef FADMM
    #include "solver_fadmm.c"
    #include "solve_KKT_fadmm.c"
    #include "FADMM_matrices.h"
#endif   
        
// c files for CPI + header for data
#ifdef CPI
    #include "solver_cp1.c"
    #include "CPI_matrices.h"
#endif  
        
// c files for CPII + header for data
#ifdef CPII
    #include "solver_cp2.c"
    #include "CPII_matrices.h"
#endif

// common c files - NOTE: maybe I should categorize them
#include "basic_operations.c"   
#include "compute_vector_to_project.c"
#include "compute_projection.c"
#include "projections.c"        
#include "compute_over_relaxation.c"          
#include "compute_projection_conjugates.c"
#include "compute_dual_residual.c"
#include "compute_primal_residual.c"
#include "compute_dual_residual_chambolle.c"     
#include "compute_primal_residual_chambolle.c"
#include "compute_array_x.c"
#include "compute_array_x_varying_stepsize.c"        
#include "compute_array_xbar.c"        
#include "compute_array_p.c"
#include "compute_array_nu.c"         
        
// main mex routine
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{ 
    
    int numiter;
    double runtime;
    double *solution;
    double *qqq;
    
    // set the solution
    plhs[0] = mxCreateDoubleMatrix(X_LEN, 1, mxREAL);
    solution = mxGetPr(plhs[0]);
    
    // set the solution
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    qqq = mxGetPr(plhs[1]);
    
    /* Do the evaluation in a subroutine */
    #ifdef ADMM
    numiter = solver_admm(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                          ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                          L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                          lam, t, y, prev_y, lx, const_vec, l, solution, q,
                          prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                          MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif
                    
    #ifdef AMA          
    numiter = solver_ama(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                         ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                         L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                         lam, t, y, prev_y, lx, const_vec, l, solution, q,
                         prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif
            
    #ifdef FAMA       
    numiter = solver_fama(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                          ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                          L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                          lam, prev_lam, lam_hat, t, y, prev_y, lx, const_vec, l, solution, q,
                          prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                          MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, ALPHA, &runtime);
    #endif
                    
    #ifdef FADMM                  
    numiter = solver_fadmm(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS, 
                           ALPHAL_ROWS, ALPHAL_COLS,  ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                           L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI, 
                           lam, prev_lam, lam_hat, t, y, prev_y, y_hat, lx, const_vec, l, solution, q,
                           prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                           MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif   
            
    #ifdef CPI
    numiter = solver_cp1(NPROX, L_ROWS, L_COLS, LI_ROWS, LI_COLS, A_ROWS, A_COLS, AI_ROWS, AI_COLS, MINVERSE_ROWS, MINVERSE_COLS, KT_ROWS, KT_COLS, KTI_ROWS, KTI_COLS, 
                         SD1_ROWS, SD1_COLS, SD1I_ROWS, SD1I_COLS, SD2_ROWS, SD2_COLS, SD2I_ROWS, SD2I_COLS, K_ROWS, K_COLS, KI_ROWS, KI_COLS, 
                         SP_ROWS, SP_COLS, SPI_ROWS, SPI_COLS, P_LEN, NU_LEN, X_LEN, P_NU_LEN, L_LEN, B_LEN, F_LEN,
                         p_prev, nu_prev, x_prev, p_nu_prev, p, nu, solution, p_nu, l, xbar, xbar_prev, b, f, L, LI, A, AI, mInverse, K, KI, KT, KTI,
                         SD1, SD1I, SD2, SD2I, SP, SPI, prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         TAU, SIGMA, THETA, PRIMALTOL, DUALTOL, MAXITER, &runtime);
    #endif            
            
    #ifdef CPII            
    numiter = solver_cp2(NPROX, L_ROWS, L_COLS, LI_ROWS, LI_COLS, A_ROWS, A_COLS, AI_ROWS, AI_COLS, KT_ROWS, KT_COLS, KTI_ROWS, KTI_COLS, 
                         SD1_ROWS, SD1_COLS, SD1I_ROWS, SD1I_COLS, SD2_ROWS, SD2_COLS, SD2I_ROWS, SD2I_COLS, K_ROWS, K_COLS, KI_ROWS, KI_COLS,
                         SP_ROWS, SP_COLS, SPI_ROWS, SPI_COLS, P_LEN, NU_LEN, X_LEN, P_NU_LEN, L_LEN, B_LEN, F_LEN, DIAGQ_LEN, p_prev, nu_prev, 
                         x_prev, p_nu_prev, diagQ, p, nu, solution, p_nu, l, xbar, xbar_prev, b, f, L, LI, A, AI, K, KI, KT, KTI, SD1, SD1I, SD2, SD2I, SP, SPI, 
                         prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         TAU, SIGMA, THETA, GAMMA, PRIMALTOL, DUALTOL, MAXITER, &runtime);
    #endif
            
    solution = x;
    qqq[0] = numiter;
         
    return;
}