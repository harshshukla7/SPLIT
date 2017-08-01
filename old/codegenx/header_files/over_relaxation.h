#ifndef _COMPUTE_OVER_RELAXATION__H__
#define _COMPUTE_OVER_RELAXATION__H__

/************************************************************************
 * function name   : compute_over_relaxation                            *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   lr : number of rows in matrix L                    *
 *                   lc : number of columns in matrix L                 *
 *                   lir : number of rows in matrix LI                  *                         
 *                   lic : number of columns in matrix LI               *
 *                   p_len : length of array arrayP                     *
 *                   l_len : length of array arrayL                     *
 *                   x_len : lenght of array arrayXbar                  *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixL, matrixLI                       *
 *                   vectors  : arrayL, arrayX, arrayP                  *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   sigma : stepsize for chambolle                     *
 ************************************************************************/


void compute_array_p(int alr,   int alc,  double matrixAlphaL[alr][alc],
                     int alir,  int alic, int matrixAlphaLI[alir][alic], 
                     int l_len, double arrayL[l_len],
                     int y_len, double arrayY[y_len], 
                     int x_len, double arrayX[x_len], double alpha);

#endif