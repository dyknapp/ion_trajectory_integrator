#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 20){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "20 inputs required.");
      }
      // C     Subroutine for Verlet integration of ion trajectory.
      //       subroutine integrate_trajectory(xx, yy, zz, vxx, vyy, vzz, : 6
      //                                       0   1   2   3    4    5
      //      &                        potential_maps, voltages, step_times_in, : 3
      //                               6               7         8
      //      &                        time_steps, dimensions, is_electrode, : 3
      //                               9           10          11
      //      &                        n_electrodes, m, q, din, maxdist, maxt, : 6 (sum: 18)
      //                               12            13 14 15   16       17
      //      &                        x_traj, y_traj, z_traj, ts, exs, eys, ezs, its) : 8
      //                               0       1       2       3   4    5    6    8
      //      & bind(c, name = "integrate_trajectory")


      // INPUTS
      const int    pts = mxGetScalar(prhs[0]);
      const int    n_particles = mxGetScalar(prhs[1]);
      mxDouble *xx  = mxGetDoubles(prhs[2]);
      mxDouble *yy  = mxGetDoubles(prhs[3]);
      mxDouble *zz  = mxGetDoubles(prhs[4]);
      mxDouble *vxx = mxGetDoubles(prhs[5]);
      mxDouble *vyy = mxGetDoubles(prhs[6]);
      mxDouble *vzz = mxGetDoubles(prhs[7]); // 8

      mxDouble *potential_maps = mxGetDoubles(prhs[8]);
      mxDouble *voltages = mxGetDoubles(prhs[9]);
      mxDouble *step_times_in = mxGetDoubles(prhs[10]); //3

      const int time_steps = mxGetScalar(prhs[11]);
      mxInt32 *dimensions = mxGetInt32s(prhs[12]);
      mxInt32 *is_electrode = mxGetInt32s(prhs[13]); // 3

      const int n_electrodes = mxGetScalar(prhs[14]);
      const double m = mxGetScalar(prhs[15]);
      const double q = mxGetScalar(prhs[16]);
      const double din = mxGetScalar(prhs[17]);
      const double maxdist = mxGetScalar(prhs[18]);
      const double maxt = mxGetScalar(prhs[19]);  // 6 (sum: 18)


      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(n_particles, pts, mxREAL);
      double *x_trajs = mxGetPr(plhs[0]);
      plhs[1] = mxCreateDoubleMatrix(n_particles, pts, mxREAL);
      double *y_trajs = mxGetPr(plhs[1]);
      plhs[2] = mxCreateDoubleMatrix(n_particles, pts, mxREAL);
      double *z_trajs = mxGetPr(plhs[2]);
      plhs[3] = mxCreateDoubleMatrix(n_particles, pts, mxREAL);
      double *tss     = mxGetPr(plhs[3]);
      plhs[4] = mxCreateDoubleMatrix(n_particles, 1,   mxREAL);
      double *itss    = mxGetPr(plhs[4]);

      mexPrintf("Calling FORTRAN code.\n");
      fly_ensemble(&pts, &n_particles, xx, yy, zz, vxx, vyy, vzz, \
                   potential_maps, voltages, step_times_in, \
                   &time_steps, dimensions, is_electrode, \
                   &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                   x_trajs, y_trajs, z_trajs, tss, itss);
      mexPrintf("FORTRAN exited successfully.\n");
}
