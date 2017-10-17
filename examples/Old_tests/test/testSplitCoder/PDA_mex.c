/*
 * Mex file to interface to generated SPLIT ADMM
 */

#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <string.h>
#include "PDA.h"
#include "splitTimer.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of input and output arguments */
    if (nrhs > 3 || nrhs < 1) 
        mexErrMsgIdAndTxt( "MATLAB:admm_mex:invalidNumInputs",
                "Usage : sol = PDA_mex(par, opt, verbose)");
    if(nrhs > 1 && !mxIsStruct(prhs[1]))
        mexErrMsgIdAndTxt( "MATLAB:admm_mex:inputNotStruct",
                "Second input must be a structure.");
    if(nlhs > 1)
        mexErrMsgIdAndTxt( "MATLAB:admm_mex:maxlhs",
                "Usage : sol = PDA_mex(par, opt, verbose)");

   // Get options if specified
    Opt opt = {1e-4, 1e-4, 1000, 10};
    int number_of_solves = 1;
    if (nrhs > 1) {
      if (mxGetField(prhs[1], 0, "primalTol")         != NULL) 
            opt.primalTol         = mxGetScalar(mxGetField(prhs[1], 0, "primalTol"));
      if (mxGetField(prhs[1], 0, "dualTol")           != NULL) 
            opt.dualTol           = mxGetScalar(mxGetField(prhs[1], 0, "dualTol"));
      if (mxGetField(prhs[1], 0, "MAXITR")            != NULL) 
            opt.MAXITR            = mxGetScalar(mxGetField(prhs[1], 0, "MAXITR"));
      if (mxGetField(prhs[1], 0, "ITR_PER_CONV_TEST") != NULL) 
            opt.ITR_PER_CONV_TEST = mxGetScalar(mxGetField(prhs[1], 0, "ITR_PER_CONV_TEST"));
      if (mxGetField(prhs[1], 0, "number_of_solves") != NULL) 
            number_of_solves = mxGetScalar(mxGetField(prhs[1], 0, "number_of_solves"));

    }

    if (nrhs > 2) { // Verbose mode
      printf("Solving with PDA using parameters:\n");
      printf("  primalTol         = %e\n", opt.primalTol);
      printf("  dualTol           = %e\n", opt.dualTol);
      printf("  MAXITR            = %i\n", opt.MAXITR);
      printf("  ITR_PER_CONV_TEST = %i\n", opt.ITR_PER_CONV_TEST);
      printf("  number_of_solves  = %i\n", number_of_solves);
    }

  // Get the input parameter
//   if (mxGetN(prhs[0])*mxGetM(prhs[0]) != par_ex_rows*par_ex_cols)
//         mexErrMsgIdAndTxt( "MATLAB:admm_mex:invalidInputSize",
//                 "Size of parameter vector incorrect");

    
    double *par = mxGetPr(prhs[0]);

  // Load data from file
  loadData();

  Sol sol;
  split_tic();
  for (int i=0; i<number_of_solves; i++) {
    // Initialize values of all variables to zero on first call
    initialize();
    solve(&sol, par, &opt);
  }
  double solve_time = split_toc() / (double)number_of_solves;

  // Fill in the matlab return structure
  const char *fnames[] = {"primal", "dual", "itr", "rDual", "rPrimal", "solve_time_ns"};
  int nFields = 6;
  plhs[0] = mxCreateStructMatrix(1, 1, nFields, fnames);

  mxArray *primal = mxCreateDoubleMatrix(nPrimal, 1, mxREAL);  
  memcpy(mxGetPr(primal), sol.primal, sizeof(double)*nPrimal);
  mxSetField(plhs[0], 0, "primal", primal);

  mxArray *dual   = mxCreateDoubleMatrix(nDual, 1, mxREAL);  
  memcpy(mxGetPr(dual), sol.dual, sizeof(double)*nDual);
  mxSetField(plhs[0], 0, "dual", dual);

  mxSetField(plhs[0], 0, "itr",     mxCreateDoubleScalar((double)(sol.itr)));
  mxSetField(plhs[0], 0, "rDual",   mxCreateDoubleScalar((double)(sol.rDual)));
  mxSetField(plhs[0], 0, "rPrimal", mxCreateDoubleScalar((double)(sol.rPrimal)));
  mxSetField(plhs[0], 0, "solve_time_ns", mxCreateDoubleScalar(solve_time));
}