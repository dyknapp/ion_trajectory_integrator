#include "mex.h"
#include "matrix.h"
#include "potential_calculation.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 5

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }
     //  subroutine nn_interpolate2D(array, odim1, odim2, ndim1, ndim2,
     // &                            new_array)
     // &                            bind(c, name="nn_interpolate2D")

      // INPUTS
      mxDouble *array               = mxGetDoubles(prhs[0]);
      const int odim1               = mxGetScalar(prhs[1]);
      const int odim2               = mxGetScalar(prhs[2]);
      const int ndim1               = mxGetScalar(prhs[3]);
      const int ndim2               = mxGetScalar(prhs[4]);

      // OUTPUTS
      plhs[0]                       = mxCreateDoubleMatrix(ndim1, ndim2, mxREAL);
      double *new_array       = mxGetPr(plhs[0]);

      nn_interpolate2D(array, &odim1, &odim2, &ndim1, &ndim2, new_array);
}