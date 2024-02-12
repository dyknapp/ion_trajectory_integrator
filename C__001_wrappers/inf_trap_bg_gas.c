#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "inf_trap_bg_gas.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     //        subroutine inf_trap_bg_gas(position, velocity, m, q
     // &              omega, depth, R, max_t, max_dist, averaging_window,
     // &              pts, four_trajectory, its)
      mxDouble *position  = mxGetDoubles(prhs[0]);
      mxDouble *velocity  = mxGetDoubles(prhs[1]);
      const double m = mxGetScalar(prhs[2]);
      const double q = mxGetScalar(prhs[3]);
      const double coll_rate = mxGetScalar(prhs[4]);
      const double room_temp = mxGetScalar(prhs[5]);
      const double omega = mxGetScalar(prhs[6]);
      const double depth = mxGetScalar(prhs[7]);
      const double R = mxGetScalar(prhs[8]);
      const double max_t = mxGetScalar(prhs[9]);
      const double max_dist = mxGetScalar(prhs[10]);
      const int averaging_window = mxGetScalar(prhs[11]);
      const int pts = mxGetScalar(prhs[12]);

      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(pts, 4, mxREAL);
      double *four_trajectory = mxGetPr(plhs[0]);
      int its, rpts, collisions;

      inf_trap_bg_gas(position, velocity, &m, &q, &coll_rate, &room_temp,\
                      &omega, &depth, &R, &max_t, &max_dist, &averaging_window, \
                      &pts, four_trajectory, &its, &rpts, &collisions);


      plhs[1] = mxCreateDoubleScalar((double)its);
      plhs[2] = mxCreateDoubleScalar((double)rpts);
      plhs[3] = mxCreateDoubleScalar((double)collisions);
}
