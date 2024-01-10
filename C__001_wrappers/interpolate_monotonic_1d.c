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
      mxDouble *f             =  mxGetDoubles(prhs[0]);
      mxDouble *t             =  mxGetDoubles(prhs[1]);
      mxDouble *new_t         =  mxGetDoubles(prhs[2]);

      int  in_length  = mxGetNumberOfElements(prhs[0]);
      int out_length  = mxGetNumberOfElements(prhs[2]);


      // OUTPUTS
      plhs[0]         = mxCreateDoubleMatrix(out_length, 1, mxREAL);
      double *new_f                 = mxGetPr(plhs[0]);

      interpolate_monotonic_1d(&in_length, &out_length, t, f, new_t, new_f);
}