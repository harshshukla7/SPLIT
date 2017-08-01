#ifndef _SOLVER_CP2__H__
#define _SOLVER_CP2__H__

/************************************************************************
 * function name   : solver_cp2                                         *  
 * -------------------------------------------------------------------- *
 * return type     : number of iterations                               *  
 * -------------------------------------------------------------------- *
 * input parameters : LOT OF STUFF                                      *
 * TODO: ADD COMMENTS                                                   *
 ************************************************************************/

int solver_cp2(int nprox, int lr, int lc, int lir, int lic, int ar, int ac, int air, int aic, int ktr, int ktc, int ktir, int ktic, 
               int sd1r, int sd1c, int sd1ir, int sd1ic, int sd2r, int sd2c, int sd2ir, int sd2ic,  int kr, int kc, int kir, int kic,
               int spr, int spc, int spir, int spic, int p_len, int nu_len, int x_len, int p_nu_len, int l_len, int b_len, int f_len, int diagQ_len,
               double arrayPprev[p_len], double arrayNuPrev[p_len], double arrayXprev[x_len], double arrayPnuPrev[p_nu_len], 
               double arrayDiagQ[diagQ_len], double arrayP[p_len], double arrayNu[nu_len],  double arrayX[x_len], double arrayPnu[p_nu_len], 
               double arrayL[l_len], double arrayXbar[x_len], double arrayXbarPrev[x_len], double arrayB[b_len], double arrayF[f_len],
               double matrixL[lr][lc], int matrixLI[lir][lic], double matrixA[ar][ac], int matrixAI[air][aic], double matrixK[kr][kc], 
               int matrixKI[kir][kic],double matrixKT[ktr][ktc], int matrixKTI[ktir][ktic], double matrixSD1[sd1r][sd1c], 
               int matrixSD1I[sd1ir][sd1ic], double matrixSD2[sd2r][sd2c], int matrixSD2I[sd2ir][sd2ic], double matrixSP[spr][spc], 
               int matrixSPI[spir][spic], char *stringArrayProxTypes[nprox], char *stringArrayNormTypes[nprox], int arrayLenOfVectors[nprox], 
               int arrayNewProxBeg[nprox], int arrayNewProxEnd[nprox], double arrayConstants[nprox], double arrayWeights[nprox], 
               double tau, double sigma, double theta, double gamma, double primaltol, double dualtol, int maxiter, double *runtime) ;
                
#endif