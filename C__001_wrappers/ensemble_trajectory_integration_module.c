#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILES ALSO
#define EXPECTED_INPUTS 19

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != EXPECTED_INPUTS){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "EXPECTED_INPUTS inputs required.");
      }

      // INPUTS
      const int particles           = mxGetScalar(prhs[0]);
      mxDouble *xs                  = mxGetDoubles(prhs[1]);
      mxDouble *ys                  = mxGetDoubles(prhs[2]);
      mxDouble *zs                  = mxGetDoubles(prhs[3]);
      mxDouble *vxs                 = mxGetDoubles(prhs[4]);
      mxDouble *vys                 = mxGetDoubles(prhs[5]);
      mxDouble *vzs                 = mxGetDoubles(prhs[6]); // 6

      mxDouble *potential_maps      = mxGetDoubles(prhs[7]);
      mxDouble *voltages            = mxGetDoubles(prhs[8]);
      mxDouble *step_times_in       = mxGetDoubles(prhs[9]); //3

      const int time_steps          = mxGetScalar(prhs[10]);
      mxInt32 *dimensions           = mxGetInt32s(prhs[11]);
      mxInt32 *is_electrode         = mxGetInt32s(prhs[12]); // 3

      const int n_electrodes        = mxGetScalar(prhs[13]);
      const double m                = mxGetScalar(prhs[14]);
      const double q                = mxGetScalar(prhs[15]);
      const double din              = mxGetScalar(prhs[16]);
      const double maxdist          = mxGetScalar(prhs[17]);
      const double maxt             = mxGetScalar(prhs[18]); // 6 (sum: 18)


      // OUTPUTS
      plhs[0]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *x_trajs         = mxGetPr(plhs[0]);
      plhs[1]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *y_trajs         = mxGetPr(plhs[1]);
      plhs[2]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *z_trajs         = mxGetPr(plhs[2]);
      plhs[3]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *tss             = mxGetPr(plhs[3]);
      plhs[4]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *exss            = mxGetPr(plhs[4]);
      plhs[5]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *eyss            = mxGetPr(plhs[5]);
      plhs[6]                       = mxCreateDoubleMatrix(particles, MAX_TRAJECTORY_POINTS, mxREAL);
      double *ezss            = mxGetPr(plhs[6]);
      plhs[7]                       = mxCreateDoubleMatrix(particles, 1, mxREAL);
      double *itss            = mxGetPr(plhs[7]);

      fly_ensemble(&particles, xs, ys, zs, vxs, vys, vzs, \
                           potential_maps, voltages, step_times_in, \
                           &time_steps, dimensions, is_electrode, \
                           &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                           x_trajs, y_trajs, z_trajs, tss, exss, eyss, ezss, itss);
}
