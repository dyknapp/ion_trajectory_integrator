#include "mex.h"
#include "matrix.h"
#include "potential_calculation.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 7

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }
     // subroutine refined_laplace(guess, bc_mask, threshold, maxits, 
     // &                           d1, d2, refinements
     // &                           its, refined)
     // & bind(c, name="iterate_laplace")

      // INPUTS
      mxDouble *guess               = mxGetDoubles(prhs[0]);
      mxDouble *bc_mask             = mxGetDoubles(prhs[1]);
      const double threshold        = mxGetScalar(prhs[2]);
      const double maxits           = mxGetScalar(prhs[3]);
      const int dim1                = mxGetScalar(prhs[4]);
      const int dim2                = mxGetScalar(prhs[5]);
      const int refinements         = mxGetScalar(prhs[6]);


      // OUTPUTS
      plhs[0]                       = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
      double *refined         = mxGetPr(plhs[0]);
      double its;

      refined_laplace(guess, bc_mask, &threshold, &maxits, &dim1, &dim2, &refinements, refined);
}
