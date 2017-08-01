#include <math.h>
#include "mex.h"
#include "matrix_ops.h"
#include "splitTimer.h"

#include "admm.h"
#include "matrix_ops.h"

/* Input Arguments */

#define PAR_IN  prhs[0]

/* Output Arguments */

#define X_OUT  plhs[0]
#define rPRIMAL_OUT  plhs[1]
#define rDUAL_OUT  plhs[2]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    if (nrhs != 1) { 
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "One input arguments required."); 
    }

    double *par = mxGetPr(PAR_IN);
    Sol sol;

    // Time the call
    split_tic();
    for(int i=0; i<1000; i++)
    {
        initialize();    // Initialize all variables to zero
        solve(&sol,par); // Solve
    }
    double toc = split_toc();
    mexPrintf("Total computation time : %e sec\n", toc/1000.0);
    mexPrintf("Computation time per iteration : %e sec\n", toc/1000.0/sol.itr);

    X_OUT = mxCreateDoubleMatrix( (mwSize)nPrimal, (mwSize)1, mxREAL); 
    double *x = mxGetPr(X_OUT);
    memcpy(x, sol.primal, sizeof(double)*nPrimal);
    
    mexPrintf("rPrimal = %g rDual = %g Num iterations = %i\n", sol.rPrimal, sol.rDual, sol.itr);

    return;    
}
