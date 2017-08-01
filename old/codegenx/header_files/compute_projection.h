#ifndef _COMPUTE_PROJECTION__H__
#define _COMPUTE_PROJECTION__H__

/************************************************************************
 * function name   : compute_projection                                 *  
 * -------------------------------------------------------------------- *
 * return type     : none                                               *  
 * -------------------------------------------------------------------- *
 * input parameters: sizes of 1D and 2D arrays                          *
 * -------------------------------------------------------------------- *
 *                   q_len : length of array arrayQ                     *        
 *                   l_len : length of array arrayL                     *
 *                   -------------------------------------------------- *
 *                   names of stringarrays and arrays                   *                     
 *                   -------------------------------------------------- *
 *                   stringArrays: stringArrayProxTypes,                *
 *                                 stringArrayNormTypes                 *
 *                   vectors     : arrayLenOfVectors                    *
 *                                 arrayNewProxBeg, arrayNewProxEnd     *
 *                                 arrayWeights,    arrayConstants      *
 *                                 arrayQ,          arrayY              *
 *                                 constants,       weights             *
 *                   -------------------------------------------------- *
 *                   miscellaneous                                      *
 *                   -------------------------------------------------- *
 *                   stepsize : step length                             *  
 *                   nprox    : number of prox operators                *   
 ************************************************************************/


void compute_projection(int nprox, char *stringArrayProxTypes[nprox],
                        char *stringArrayNormTypes[nprox],
                        int arrayLenOfVectors[nprox],
                        int arrayNewProxBeg[nprox], 
                        int arrayNewProxEnd[nprox], 
                        double arrayConstants[nprox], double arrayWeights[nprox],
                        int q_len, double arrayQ[q_len],
                        int y_len, double arrayY[y_len], double stepsize);

#endif