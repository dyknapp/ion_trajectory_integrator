#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 21){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "20 inputs required.");
      }


      // INPUTS
      const double record_interval = mxGetScalar(prhs[0]);
      const int            interps = mxGetScalar(prhs[1]);
      const int        n_particles = mxGetScalar(prhs[2]);
      mxDouble                *xx  = mxGetDoubles(prhs[3]);
      mxDouble                *yy  = mxGetDoubles(prhs[4]);
      mxDouble                *zz  = mxGetDoubles(prhs[5]);
      mxDouble                *vxx = mxGetDoubles(prhs[6]);
      mxDouble                *vyy = mxGetDoubles(prhs[7]);
      mxDouble                *vzz = mxGetDoubles(prhs[8]); // 8

      mxDouble     *potential_maps = mxGetDoubles(prhs[9]);
      mxDouble           *voltages = mxGetDoubles(prhs[10]);
      mxDouble      *step_times_in = mxGetDoubles(prhs[11]); //3

      const int         time_steps = mxGetScalar(prhs[12]);
      mxInt32          *dimensions = mxGetInt32s(prhs[13]);
      mxInt32        *is_electrode = mxGetInt32s(prhs[14]); // 3

      const int       n_electrodes = mxGetScalar(prhs[15]);
      mxDouble                 *ms = mxGetDoubles(prhs[16]);
      mxDouble                 *qs = mxGetDoubles(prhs[17]);
      const double             din = mxGetScalar(prhs[18]);
      const double         maxdist = mxGetScalar(prhs[19]);
      const double            maxt = mxGetScalar(prhs[20]);  // 6 (sum: 18)


      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(n_particles, interps, mxREAL);
      double *x_trajs = mxGetPr(plhs[0]);
      plhs[1] = mxCreateDoubleMatrix(n_particles, interps, mxREAL);
      double *y_trajs = mxGetPr(plhs[1]);
      plhs[2] = mxCreateDoubleMatrix(n_particles, interps, mxREAL);
      double *z_trajs = mxGetPr(plhs[2]);
      plhs[3] = mxCreateDoubleMatrix(1, interps, mxREAL);
      double *tss     = mxGetPr(plhs[3]);

      double its;
      mexPrintf("Calling FORTRAN code.\n");
      fly_cloud(&record_interval, &interps, &n_particles,\
                xx, yy, zz, vxx, vyy, vzz, \
                potential_maps, voltages, step_times_in, \
                &time_steps, dimensions, is_electrode, \
                &n_electrodes, ms, qs, &din, &maxdist, &maxt, \
                x_trajs, y_trajs, z_trajs, tss, &its);
      mexPrintf("FORTRAN exited successfully.\n");

      plhs[4] = mxCreateDoubleScalar(its);
}
