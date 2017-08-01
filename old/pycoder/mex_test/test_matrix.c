#include <math.h>
#include "mex.h"
#include "matrix_ops.h"

/* Input Arguments */

#define X_IN  prhs[0]
#define C_IN  prhs[1]
#define FUNC  prhs[2]

/* Output Arguments */

#define YP_OUT  plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    if (nrhs != 3) { 
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "Three input arguments required."); 
    }

    int m = mxGetM(X_IN); 
    int n = mxGetN(X_IN);
    YP_OUT = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n, mxREAL); 
    double *yp = mxGetPr(YP_OUT);
    
    double c   = mxGetScalar(C_IN); 
    double *y  = mxGetPr(X_IN);
    int func = mxGetScalar(FUNC);

    switch(func) {
        case 1:
            prox_norm_one(yp, y, c, m);
            break;
        case 2:
            prox_norm_two(yp, y, c, m);
            break;
        case 3:
            proj_normball_two(yp, y, c, m);
            break;
        case 4:
            proj_normball_inf(yp, y, c, m);
            break;
        case 5:
            proj_secondOrderCone(yp, y, m);
            break;
        case 6:
            proj_negative(yp, y, m);
            break;
        case 7:
            proj_positive(yp, y, m);
            break;
        }
return;    
}
