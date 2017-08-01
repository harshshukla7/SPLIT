#ifndef _SOLVER_FAMA__H__
#define _SOLVER_FAMA__H__

/************************************************************************
 * function name   : solver_fama                                        *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters : LOT OF STUFF                                      *
 * TODO: ADD COMMENTS                                                   *
 ************************************************************************/

int solver_fama(int mr, int mc, int lr, int lc, int mir, int mic, int lir, int lic, int slr, int slc, int slir, int slic, int alr, int alc, 
                int alir, int alic, int nprox, int lam_len, int t_len, int y_len, int cv_len, int x_len, int q_len, int l_len, int lx_len, 
                double matrixL[lr][lc],    int matrixLI[lr][lc], 
                double matrixM[mr][mc],    int matrixMI[mir][mic], 
                double matrixSL[slr][slc], int matrixSLI[slir][slic],
                double matrixAlphaL[alr][alc], int matrixAlphaLI[alir][alic],
                double arrayLam[lam_len],  double arrayPrevLam[lam_len], double arrayLamHat[lam_len],  
                double arrayT[t_len], double arrayY[y_len], double arrayPrevY[y_len], double arrayLx[lx_len], 
                double arrayConstVec[cv_len], double arrayL[l_len], double arrayX[x_len], double arrayQ[q_len],
                char *stringArrayProxTypes[nprox], char *stringArrayNormTypes[nprox],
                int arrayLenOfVectors[nprox], int arrayNewProxBeg[nprox], int arrayNewProxEnd[nprox], 
                double arrayConstants[nprox], double arrayWeights[nprox],
                int maxiter, int num_opt_var, double stepsize, double dualtol, double alpha, double *runtime);

#endif