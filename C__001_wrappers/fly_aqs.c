#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 24

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }

      // INPUTS
      mxDouble *amp_scales          = mxGetDoubles(prhs[0]);
      const int a_res               = mxGetScalar(prhs[1]);
      mxDouble *off_scales          = mxGetDoubles(prhs[2]);
      const int o_res               = mxGetScalar(prhs[3]);
      const int reps                = mxGetScalar(prhs[4]);
      const double xx               = mxGetScalar(prhs[5]);
      const double yy               = mxGetScalar(prhs[6]);
      const double zz               = mxGetScalar(prhs[7]);
      const double vxx              = mxGetScalar(prhs[8]);
      const double vyy              = mxGetScalar(prhs[9]);
      const double vzz              = mxGetScalar(prhs[10]); // 11

      mxDouble *rf_potential        = mxGetDoubles(prhs[11]);
      mxDouble *endcap_potential    = mxGetDoubles(prhs[12]);
      mxDouble *voltages            = mxGetDoubles(prhs[13]);
      mxDouble *step_times_in       = mxGetDoubles(prhs[14]); // 4

      const int time_steps          = mxGetScalar(prhs[15]);
      mxInt32 *dimensions           = mxGetInt32s(prhs[16]);
      mxInt32 *is_electrode         = mxGetInt32s(prhs[17]); // 3

      const int n_electrodes        = mxGetScalar(prhs[18]);
      const double m                = mxGetScalar(prhs[19]);
      const double q                = mxGetScalar(prhs[20]);
      const double din              = mxGetScalar(prhs[21]);
      const double maxdist          = mxGetScalar(prhs[22]);
      const double maxt             = mxGetScalar(prhs[23]); // 6 (sum: 18)


      // OUTPUTS
      plhs[0]                       = mxCreateDoubleMatrix(a_res, o_res, mxREAL);
      double *lifetimes       = mxGetPr(plhs[0]);

      fly_aqs(amp_scales, &a_res, off_scales, &o_res, &reps, \
                              &xx, &yy, &zz, &vxx, &vyy, &vzz, \
                              rf_potential, endcap_potential, \
                              voltages, step_times_in, \
                              &time_steps, dimensions, is_electrode, \
                              &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                              lifetimes);
}
