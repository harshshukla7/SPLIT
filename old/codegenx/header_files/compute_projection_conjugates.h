#ifndef _COMPUTE_PROJECTION_CONJUGATES__H__
#define _COMPUTE_PROJECTION_CONJUGATES__H__

/************************************************************************
 * function name   : compute_projection_conjugates                      *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   p_len : length of array arrayP                     *        
 *                   -------------------------------------------------- *
 *                   names of stringarrays and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   stringArrays: stringArrayProxTypes,                *
 *                                 stringArrayNormTypes                 *
 *                   vectors     : arrayLenOfVectors, arrayP            *
 *                                 arrayNewProxBeg, arrayNewProxEnd     *
 *                                 arrayWeights,    arrayConstants      *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   stepsize : step length: reduntant here: modify it  *  
 *                   nprox    : number of prox operators                *   
 ************************************************************************/


void compute_projection_conjugates(int nprox, char *stringArrayProxTypes[nprox],
                                   char *stringArrayNormTypes[nprox],
                                   int arrayLenOfVectors[nprox],
                                   int arrayNewProxBeg[nprox], 
                                   int arrayNewProxEnd[nprox], 
                                   double arrayConstants[nprox], double arrayWeights[nprox],
                                   int p_len, double arrayP[p_len]);

#endif