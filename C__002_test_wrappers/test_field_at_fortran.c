#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"
#define MAX_TRAJECTORY_POINTS 131072
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 11){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "18 inputs required.");
      }


      // INPUTS
      const double t = mxGetScalar(prhs[0]);
      mxDouble *voltages = mxGetDoubles(prhs[1]);
      mxDouble *step_times = mxGetDoubles(prhs[2]);
      const int n_electrodes = mxGetScalar(prhs[3]);
      const int time_steps = mxGetScalar(prhs[4]);
      const double xx = mxGetScalar(prhs[5]);
      const double yy = mxGetScalar(prhs[6]);
      const double zz = mxGetScalar(prhs[7]);
      const double din = mxGetScalar(prhs[8]);
      mxDouble *potential_maps = mxGetDoubles(prhs[9]);
      mxInt32 *dimensions = mxGetInt32s(prhs[10]);

      // OUTPUTS
      double ex, ey, ez;
      field_at(&t, voltages, step_times, \
           &n_electrodes, &time_steps,\
           &xx, &yy, &zz, &din, potential_maps, dimensions,\
           &ex, &ey, &ez);

      plhs[0] = mxCreateDoubleScalar(ex);
      plhs[1] = mxCreateDoubleScalar(ey);
      plhs[2] = mxCreateDoubleScalar(ez);
}