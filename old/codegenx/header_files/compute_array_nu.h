#ifndef _COMPUTE_ARRAY_NU__H__
#define _COMPUTE_ARRAY_NU__H__

/************************************************************************
 * function name   : compute_array_nu                                   *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   ar : number of rows in matrix A                    *
 *                   ac : number of columns in matrix A                 *
 *                   air : number of rows in matrix AI                  *                         
 *                   aic : number of columns in matrix AI               *
 *                   nu_len : length of array arrayNu                   *
 *                   b_len : length of array arrayB                     *
 *                   x_len : lenght of array arrayXbar                  *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixA, matrixAI                       *
 *                   vectors  : arrayB, arrayXbar, arrayNu              *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   sigma : stepsize for chambolle                     *
 ************************************************************************/


void compute_array_nu(int ar, int ac, double matrixA[ar][ac],
                      int air, int aic, int matrixAI[air][aic], 
                      int b_len, double arrayB[b_len],
                      int nu_len, double arrayNu[nu_len], 
                      int x_len, double arrayXbar[x_len], double sigma);

#endif