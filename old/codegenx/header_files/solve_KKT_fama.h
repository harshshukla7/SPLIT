#ifndef _SOLVE_KKT_FAMA__H__
#define _SOLVE_KKT_FAMA__H__

/************************************************************************
 * function name   : solve_KKT_fama                                     *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   mr : number of rows in matrix M                    *
 *                   mc : number of columns in matrix M                 *
 *                   mir : number of rows in matrix MI                  *                         
 *                   mic : number of columns in matrix MI               *
 *                   lam_len : length of array arrayLam                 *
 *                   t_len: length of array arrayT                      *
 *                   y_len : length of array arrayY                     *        
 *                   cv_len : length of array arrayConstVec             *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixM, matrixMI                       *
 *                   vectors  : arrayLam, arrayY, arrayConstVec         *
 ************************************************************************/

void solve_KKT_fama(int mr, int mc, double matrixM[mr][mc],
                    int mir, int mic, int matrixMI[mir][mic],
                    int t_len, double arrayT[t_len],   
                    int lam_len, double arrayLamHat[lam_len], 
                    int y_len, double arrayY[y_len],
                    int cv_len, double arrayConstVec[cv_len]);

#endif