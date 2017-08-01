#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// define macros for solvers
// #define ADMM
// #define AMA 
// #define FAMA
// #define FADMM
// #define CPI
#define CPII

// include header files: solver + data
#ifdef ADMM
    #include "../header_files/solver_admm.h"    
    #include "ADMM_matrices.h"
#endif    
        
#ifdef AMA
    #include "../header_files/solver_ama.h" 
    #include "AMA_matrices.h"
#endif
        
#ifdef FAMA
    #include "../header_files/solver_fama.h"     
    #include "FAMA_matrices.h"
#endif
        
#ifdef FADMM
    #include "../header_files/solver_fadmm.h" 
    #include "FADMM_matrices.h"
#endif
        
#ifdef CPI
    #include "../header_files/solver_cp1.h"     
    #include "CPI_matrices.h"
#endif        

#ifdef CPII
    #include "../header_files/solver_cp2.h"     
    #include "CPII_matrices.h"
#endif
        
// common header file
#include "../header_files/basic_operations.h"

int main(void)
{
    double runtime;
    int numiter;
    
    #ifdef ADMM
    numiter = solver_admm(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                          ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                          L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                          lam, t, y, prev_y, lx, const_vec, l, x, q,
                          prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                          MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif
                    
    #ifdef AMA          
    numiter = solver_ama(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                         ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                         L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                         lam, t, y, prev_y, lx, const_vec, l, x, q,
                         prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif
            
    #ifdef FAMA       
    numiter = solver_fama(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS,
                          ALPHAL_ROWS, ALPHAL_COLS, ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                          L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI,
                          lam, prev_lam, lam_hat, t, y, prev_y, lx, const_vec, l, x, q,
                          prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                          MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, ALPHA, &runtime);
    #endif
                    
    #ifdef FADMM                  
    numiter = solver_fadmm(MATRIX_ROWS, MATRIX_COLS, L_ROWS, L_COLS, MATRIXI_ROWS, MATRIXI_COLS, LI_ROWS, LI_COLS, SCALEDL_ROWS, SCALEDL_COLS, SCALEDLI_ROWS, SCALEDLI_COLS, 
                           ALPHAL_ROWS, ALPHAL_COLS,  ALPHALI_ROWS, ALPHALI_COLS, NPROX, LAM_LEN, T_LEN, Y_LEN, CONST_VEC_LEN, X_LEN, Q_LEN, L_LEN, LX_LEN,
                           L, LI, matrix, matrixI, scaledL, scaledLI, alphaL, alphaLI, 
                           lam, prev_lam, lam_hat, t, y, prev_y, y_hat, lx, const_vec, l, x, q,
                           prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                           MAXITER, NUM_OPT_VAR, STEPSIZE, DUALTOL, PRIMALTOL, ALPHA, &runtime);
    #endif   
                   
    #ifdef CPI
    numiter = solver_cp1(NPROX, L_ROWS, L_COLS, LI_ROWS, LI_COLS, A_ROWS, A_COLS, AI_ROWS, AI_COLS, MINVERSE_ROWS, MINVERSE_COLS, KT_ROWS, KT_COLS, KTI_ROWS, KTI_COLS, 
                         SD1_ROWS, SD1_COLS, SD1I_ROWS, SD1I_COLS, SD2_ROWS, SD2_COLS, SD2I_ROWS, SD2I_COLS, K_ROWS, K_COLS, KI_ROWS, KI_COLS, 
                         SP_ROWS, SP_COLS, SPI_ROWS, SPI_COLS, P_LEN, NU_LEN, X_LEN, P_NU_LEN, L_LEN, B_LEN, F_LEN,
                         p_prev, nu_prev, x_prev, p_nu_prev, p, nu, x, p_nu, l, xbar, xbar_prev, b, f, L, LI, A, AI, mInverse, K, KI, KT, KTI,
                         SD1, SD1I, SD2, SD2I, SP, SPI, prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         TAU, SIGMA, THETA, PRIMALTOL, DUALTOL, MAXITER, &runtime);
    #endif            
            
    #ifdef CPII            
    numiter = solver_cp2(NPROX, L_ROWS, L_COLS, LI_ROWS, LI_COLS, A_ROWS, A_COLS, AI_ROWS, AI_COLS, KT_ROWS, KT_COLS, KTI_ROWS, KTI_COLS, 
                         SD1_ROWS, SD1_COLS, SD1I_ROWS, SD1I_COLS, SD2_ROWS, SD2_COLS, SD2I_ROWS, SD2I_COLS, K_ROWS, K_COLS, KI_ROWS, KI_COLS,
                         SP_ROWS, SP_COLS, SPI_ROWS, SPI_COLS, P_LEN, NU_LEN, X_LEN, P_NU_LEN, L_LEN, B_LEN, F_LEN, DIAGQ_LEN, p_prev, nu_prev, 
                         x_prev, p_nu_prev, diagQ, p, nu, x, p_nu, l, xbar, xbar_prev, b, f, L, LI, A, AI, K, KI, KT, KTI, SD1, SD1I, SD2, SD2I, SP, SPI, 
                         prox_types, norm_types, len_of_vectors, new_prox_beg, new_prox_end, constants, weights,
                         TAU, SIGMA, THETA, GAMMA, PRIMALTOL, DUALTOL, MAXITER, &runtime);
    #endif  
    
    // print number of iterations and runtime to screen
    //printf("\n------------ STATS -----------\n");
    printf("Number of iterations: %d\n", numiter);            
    printf("Elapsed time is: %f ms\n\n", runtime);
    
    // print solution to screen
    printf("----------- SOLUTION ----------\n");
    print_the_vector(x, X_LEN);
            
    return 0;
}