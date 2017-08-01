#ifndef _COMPUTE_DUAL_RESIDUAL_CHAMBOLLE_H__
#define _COMPUTE_DUAL_RESIDUAL_CHAMBOLLE_H__

/************************************************************************
 * function name   : compute_dual_residual_chambolle                    *  
 * -------------------------------------------------------------------- *
 * return type     : double: value of dual residual                      *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   sr : number of rows in matrix SCALEDK              *
 *                   sc : number of columns in matrix SCALEDK           *
 *                   sir : number of rows in matrix SCALEDKI            *
 *                   sic : number of columns in matrix SCALEDKI         *
 *                   x_len : length of array arrayX                     *
 *                   p_nu_len : length of array arrayPnu                *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices: matrixSK,  matrixSKI                     *
 *                   arrays  : arrayPnu,  arrayPnuPrev,                 *
 *                             arrayXbar, arrayXPrev                    *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   sigma : stepsize for chambolle                     *
 ************************************************************************/

double compute_dual_residual_chambolle(int kr,    int kc,    double matrixK[kr][kc],
                                       int kir,   int kic,   int matrixKI[kir][kic],
                                       int sd1r,  int sd1c,  double matrixSD1[sd1r][sd1c],
                                       int sd1ir, int sd1ic, int matrixSD1I[sd1ir][sd1ic],
                                       int sd2r,  int sd2c,  double matrixSD2[sd2r][sd2c],
                                       int sd2ir, int sd2ic, int matrixSD2I[sd2ir][sd2ic], 
                                       int x_len, double arrayX[x_len], double arrayXbarPrev[x_len], double arrayXprev[x_len],
                                       int p_nu_len, double arrayPnu[p_nu_len], double arrayPnuPrev[p_nu_len], double sigma);

#endif