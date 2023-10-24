#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 19){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "19 inputs required.");
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
      const double xx = mxGetScalar(prhs[0]);
      const double yy = mxGetScalar(prhs[1]);
      const double zz = mxGetScalar(prhs[2]);
      const double vxx = mxGetScalar(prhs[3]);
      const double vyy = mxGetScalar(prhs[4]);
      const double vzz = mxGetScalar(prhs[5]); // 6

      mxDouble *potential_maps = mxGetDoubles(prhs[6]);
      mxDouble *voltages = mxGetDoubles(prhs[7]);
      mxDouble *step_times_in = mxGetDoubles(prhs[8]); //3

      const int time_steps = mxGetScalar(prhs[9]);
      mxInt32 *dimensions = mxGetInt32s(prhs[10]);
      mxInt32 *is_electrode = mxGetInt32s(prhs[11]); // 3

      const int n_electrodes = mxGetScalar(prhs[12]);
      const double m = mxGetScalar(prhs[13]);
      const double q = mxGetScalar(prhs[14]);
      const double din = mxGetScalar(prhs[15]);
      const double maxdist = mxGetScalar(prhs[16]);
      const double maxt = mxGetScalar(prhs[17]);  // 6 (sum: 18)
      const double data_t_interval = mxGetScalar(prhs[18]);

      size_t out_size = ceil(maxt / data_t_interval);
      int out_length = (int)out_size;


      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(out_size, 1, mxREAL);
      double *x_traj = mxGetPr(plhs[0]);
      plhs[1] = mxCreateDoubleMatrix(out_size, 1, mxREAL);
      double *y_traj = mxGetPr(plhs[1]);
      plhs[2] = mxCreateDoubleMatrix(out_size, 1, mxREAL);
      double *z_traj = mxGetPr(plhs[2]);
      plhs[3] = mxCreateDoubleMatrix(out_size, 1, mxREAL);
      double *ts = mxGetPr(plhs[3]);
      double its;
      double recorded;

      mexPrintf("Calling FORTRAN code.\n");
      integrate_trajectory_lite(&xx, &yy, &zz, &vxx, &vyy, &vzz, \
                           potential_maps, voltages, step_times_in, \
                           &time_steps, dimensions, is_electrode, \
                           &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                           &data_t_interval, &out_length,
                           x_traj, y_traj, z_traj, ts, &its, &recorded);
      mexPrintf("FORTRAN exited successfully.\n");
      plhs[4] = mxCreateDoubleScalar(its);
      plhs[5] = mxCreateDoubleScalar(recorded);
}
