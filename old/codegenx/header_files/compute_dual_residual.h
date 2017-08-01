#ifndef _COMPUTE_DUAL_RESIDUAL__H__
#define _COMPUTE_DUAL_RESIDUAL__H__

/************************************************************************
 * function name   : compute_dual_residual                              *  
 * -------------------------------------------------------------------- *
 * return type     : valueof dual residual: rDual                       *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   lr : number of rows in matrix L                    *
 *                   lc : number of columns in matrix L                 *
 *                   lir : number of rows in matrix LI                  *                         
 *                   lic : number of columns in matrix LI               *
 *                   q_len : length of array arrayQ                     *        
 *                   y_len : length of array arrayY                     *
 *                   x_len : length of array arrayX                     *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixL, matrixLI                       *
 *                   vectors  : arrayY, arrayPrevY, arrayL, arrayX      *
 ************************************************************************/

double compute_dual_residual(int lr, int lc, double matrixL[lr][lc],
                             int lir, int lic, int matrixLI[lir][lic], 
                             int y_len, double arrayY[y_len],
                             int l_len, double arrayL[l_len],
                             int x_len, double arrayX[x_len]);

#endif