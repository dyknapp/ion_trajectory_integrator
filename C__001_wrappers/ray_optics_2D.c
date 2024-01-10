#include "mex.h"
#include "matrix.h"
#include "ion_optics.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      if(nrhs != 13){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "13 inputs required.");
      }

     //  subroutine ray_optics_2D(xx, yy, vxx, vyy,
     //                             0   1    2    3
     // &                        potential_maps, voltages, dimensions,
     //                                        4         5           6
     // &                        n_electrodes, m, q, din, maxdist, maxt,
     //                                      7  8  9   10       11    12
     // &                        x_traj, y_traj, ts, its)

      // INPUTS
      const double xx          =  mxGetScalar(prhs[0] );
      const double yy          =  mxGetScalar(prhs[1] );
      const double vxx         =  mxGetScalar(prhs[2] );
      const double vyy         =  mxGetScalar(prhs[3] );
      mxDouble *potential_maps = mxGetDoubles(prhs[4] );
      mxDouble *voltages       = mxGetDoubles(prhs[5] );
      mxInt32 *dimensions      =  mxGetInt32s(prhs[6] );
      const int n_electrodes   =  mxGetScalar(prhs[7] );
      const double m           =  mxGetScalar(prhs[8] );
      const double q           =  mxGetScalar(prhs[9] );
      const double din         =  mxGetScalar(prhs[10]);
      const double maxdist     =  mxGetScalar(prhs[11]);
      const double maxt        =  mxGetScalar(prhs[12]);


      // OUTPUTS
      plhs[0]        = mxCreateDoubleMatrix(MAX_TRAJECTORY_POINTS, 1, mxREAL);
      double *x_traj = mxGetPr(plhs[0]);
      plhs[1]        = mxCreateDoubleMatrix(MAX_TRAJECTORY_POINTS, 1, mxREAL);
      double *y_traj = mxGetPr(plhs[1]);
      plhs[2]        = mxCreateDoubleMatrix(MAX_TRAJECTORY_POINTS, 1, mxREAL);
      double *ts     = mxGetPr(plhs[2]);
      double its;

      mexPrintf("Calling FORTRAN code.\n");
      ray_optics_2D(&xx, &yy, &vxx, &vyy, \
                    potential_maps, voltages,\
                    dimensions, &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                    x_traj, y_traj, ts, &its);
      mexPrintf("FORTRAN exited successfully.\n");
      plhs[3] = mxCreateDoubleScalar(its);
}
