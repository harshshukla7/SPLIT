#ifndef _COMPUTE_PRIMAL_RESIDUAL__H__
#define _COMPUTE_PRIMAL_RESIDUAL__H__

/************************************************************************
 * function name   : compute_primal_residual                            *  
 * -------------------------------------------------------------------- *
 * return type     : value of primal residual: rPrimal                  *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   lr : number of rows in matrix L                    *
 *                   lc : number of columns in matrix L                 *
 *                   lir : number of rows in matrix LI                  *                         
 *                   lic : number of columns in matrix LI               *
 *                   y_len : length of array arrayY                     *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   --------------------------------------------!------ *
 *                   matrices : matrixL, matrixLI                       *
 *                   vectors  : arrayY, arrayPrevY, arrayL              *
 ************************************************************************/


double compute_primal_residual(int slr, int slc, double matrixSL[slr][slc],
                               int slir, int slic, int matrixSLI[slir][slic], 
                               int y_len, double arrayY[y_len], double arrayPrevY[y_len],
                               int l_len, double arrayL[l_len]);

#endif