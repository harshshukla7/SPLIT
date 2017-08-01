#ifndef _COMPUTE_PRIMAL_RESIDUAL_CHAMBOLLE_H__
#define _COMPUTE_PRIMAL_RESIDUAL_CHAMBOLLE_H__

/************************************************************************
 * function name   : compute_primal_residual_chambolle                  *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   x_len : length of array arrayX                     *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   arrays  : arrayX, arrayPrevX                       *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   tautheta : stepsize for chambolle: 1/tau*theta     *
 ************************************************************************/

double compute_primal_residual_chambolle(int spr,   int spc,  double matrixSP[spr][spc],
                                         int spir,  int spic, int matrixSPI[spir][spic],
                                         int x_len, double arrayX[x_len], double arrayPrevX[x_len], double tau);

#endif