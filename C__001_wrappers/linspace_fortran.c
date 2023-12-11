#include "mex.h"
#include "matrix.h"
#include "utils.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 3){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "2 inputs required.");
      }
      // INPUTS
      double    a             =  mxGetScalar(prhs[0]);
      double    b             =  mxGetScalar(prhs[1]);
      int       n             =  mxGetScalar(prhs[2]);


      // OUTPUTS
      plhs[0]         = mxCreateDoubleMatrix(n, 1, mxREAL);
      double *result          = mxGetPr(plhs[0]);

      linspace_fortran(&a, &b, &n, result);
}