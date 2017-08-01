#ifndef _COMPUTE_ARRAY_P__H__
#define _COMPUTE_ARRAY_P__H__

/************************************************************************
 * function name   : compute_array_p                                    *  
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


void compute_array_p(int lr, int lc, double matrixL[lr][lc],
                     int lir, int lic, int matrixLI[lir][lic], 
                     int l_len, double arrayL[l_len],
                     int p_len, double arrayP[p_len], 
                     int x_len, double arrayXbar[x_len], double sigma);

#endif