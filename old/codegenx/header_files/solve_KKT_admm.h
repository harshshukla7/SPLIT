#ifndef _SOLVE_KKT_ADMM__H__
#define _SOLVE_KKT_ADMM__H__

/************************************************************************
 * function name   : solve_KKT_admm                                     *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   mr : number of rows in matrix M                    *
 *                   mc : number of columns in matrix M                 *
 *                   mir : number of rows in matrix MI                  *                         
 *                   mic : number of columns in matrix MI               *
 *                   t_len: length of array arrayT                      *
 *                   lam_len : length of array arrayLam                 *
 *                   y_len : length of array arrayY                     *        
 *                   cv_len : length of array arrayConstVec             *
 *                   -------------------------------------------------- *
 *                   names of the matrices and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   matrices : matrixM, matrixMI                       *
 *                   vectors  : arrayLam, arrayY,                       *
 *                              arrayT, arrayConstVec                   *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   stepsize : step length                             *
 ************************************************************************/


void solve_KKT_admm(int mr, int mc, double matrixM[mr][mc],
                    int mir, int mic, int matrixMI[mir][mic], 
                    int lam_len, double arrayLam[lam_len],
                    int t_len, double arrayT[t_len], 
                    int y_len, double arrayY[y_len],
                    int cv_len, double arrayConstVec[cv_len], double stepsize);

#endif