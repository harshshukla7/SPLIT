#ifndef _COMPUTE_VECTOR_TO_PROJECT__H__
#define _COMPUTE_VECTOR_TO_PROJECT__H__

/************************************************************************
 * function name   : compute_vector_to_project                          *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   lx_len  : length of array arrayLx                  *
 *                   lam_len : length of array arrayLam                 *
 *                   q_len   : length of array arrayQ                   *        
 *                   l_len   : length of array arrayL                   *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixL, matrixLI                       *
 *                   vectors  : arrayLam, arrayQ, arrayL                *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   stepsize : step length                             *
 ************************************************************************/


void compute_vector_to_project(int lam_len, double arrayLam[lam_len], 
                               int q_len,   double arrayQ[q_len],
                               int l_len,   double arrayL[l_len], 
                               int lx_len,  double arrayLx[lx_len],
                               int x_len,   double arrayX[x_len], double stepsize);

#endif