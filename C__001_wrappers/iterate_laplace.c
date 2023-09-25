#include "mex.h"
#include "matrix.h"
#include "potential_calculation.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 6

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }
     //  subroutine iterate_laplace(guess, bc_mask, threshold, maxits, 
     // &                           d1,d2,
     // &                           its, refined)
     // & bind(c, name="iterate_laplace")

      // INPUTS
      mxDouble *guess               = mxGetDoubles(prhs[0]);
      mxDouble *bc_mask             = mxGetDoubles(prhs[1]);
      const double threshold        = mxGetScalar(prhs[2]);
      const double maxits           = mxGetScalar(prhs[3]);
      const int dim1                = mxGetScalar(prhs[4]);
      const int dim2                = mxGetScalar(prhs[5]);


      // OUTPUTS
      plhs[1]                       = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
      double *refined         = mxGetPr(plhs[1]);
      double its;

      iterate_laplace(guess, bc_mask, &threshold, &maxits, &dim1, &dim2, &its, refined);

      plhs[0] = mxCreateDoubleScalar(its);
}
