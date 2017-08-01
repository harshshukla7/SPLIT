#ifndef _COMPUTE_ARRAY_XBAR__H__
#define _COMPUTE_ARRAY_XBAR__H__

/************************************************************************
 * function name   : compute_array_xbar                                 *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   x_len : lenght of array arrayX                     *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   arrays  : arrayX, arrayPrevX, arrayXbar            *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   theta : stepsize for chambolle                     *
 ************************************************************************/


void compute_array_xbar(int x_len, double arrayX[x_len], double arrayXbar[x_len], 
                        double arrayPrevX[x_len], double theta);

#endif