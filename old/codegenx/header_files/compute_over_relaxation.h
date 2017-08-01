#ifndef _COMPUTE_OVER_RELAXATION__H__
#define _COMPUTE_OVER_RELAXATION__H__

/************************************************************************
 * function name   : compute_array_p                                    *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   alr  : number of rows in matrix alpha*L            *
 *                   alc  : number of columns in matrix alpha*L         *
 *                   alir : number of rows in matrix LI                 *                         
 *                   alic : number of columns in matrix LI              *
 *                   l_len : length of array arrayP                     *
 *                   y_len : length of array arrayL                     *
 *                   x_len : lenght of array arrayXbar                  *
 *                   lx_len : lenght of array arrayXbar                 *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixL, matrixLI                       *
 *                   vectors  : arrayL, arrayX, arrayP                  *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   alpha : parameter for over - relaxation            *
 ************************************************************************/


void compute_over_relaxation(int alr,    int alc,  double matrixAlphaL[alr][alc],
                             int alir,   int alic, int matrixAlphaLI[alir][alic], 
                             int lx_len, double arrayLx[lx_len],
                             int l_len,  double arrayL[l_len],
                             int y_len,  double arrayY[y_len], 
                             int x_len,  double arrayX[x_len], double alpha);

#endif