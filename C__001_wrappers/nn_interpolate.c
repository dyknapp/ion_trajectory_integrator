#include "mex.h"
#include "matrix.h"
#include "potential_calculation.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 3

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }
      // nn_interpolate(array, old_samples, new_samples, new_array)

      // INPUTS
      mxDouble *array               = mxGetDoubles(prhs[0]);
      const int old_samples         = mxGetScalar(prhs[1]);
      const int new_samples         = mxGetScalar(prhs[2]);

      // OUTPUTS
      plhs[0]                       = mxCreateDoubleMatrix(new_samples, 1, mxREAL);
      double *new_array       = mxGetPr(plhs[0]);

      nn_interpolate(array, &old_samples, &new_samples, new_array);
}