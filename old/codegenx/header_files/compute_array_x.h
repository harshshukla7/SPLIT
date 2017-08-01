#ifndef _COMPUTE_ARRAY_X__H__
#define _COMPUTE_ARRAY_X__H__

/************************************************************************
 * function name   : compute_array_x                                    *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   mir : number of rows in matrix mInverse            *
 *                   mic : number of columns in matrix mInverse         *
 *                   ktr : number of rows in matrix KT                  *                         
 *                   ktc : number of columns in matrix KT               *
 *                   p_nu_len : length of array arrayPnu                *
 *                   x_len : lenght of array arrayX                     *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixM, matrixKT, matrixKTI            *
 *                   vectors  : arrayX, arrayPnu                        *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   tau : stepsize for chambolle                       *
 *                   stepsize : redundant - modify this function        *
 ************************************************************************/


void compute_array_x(int mir,   int mic,  double matrixMI[mir][mic],
                     int ktr,   int ktc,  double matrixKT[ktr][ktc],
                     int ktir,  int ktic, int matrixKTI[ktir][ktic], 
                     int p_nu_len, double arrayPnu[p_nu_len], int x_len, double arrayX[x_len], 
                     int f_len, double arrayF[f_len], double tau);

#endif