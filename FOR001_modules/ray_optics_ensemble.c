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
      mxDouble *xxs            =  mxGetScalar(prhs[0] );
      mxDouble *yys            =  mxGetScalar(prhs[1] );
      mxDouble *vxxs           =  mxGetScalar(prhs[2] );
      mxDouble *vyys           =  mxGetScalar(prhs[3] );
      mxDouble *potential_maps = mxGetDoubles(prhs[4] );
      mxDouble *voltagess      = mxGetDoubles(prhs[5] );
      mxInt32 *dimensions      =  mxGetInt32s(prhs[6] );
      const int n_electrodes   =  mxGetScalar(prhs[7] );
      mxDouble *ms             =  mxGetScalar(prhs[8] );
      mxDouble *qs             =  mxGetScalar(prhs[9] );
      const double din         =  mxGetScalar(prhs[10]);
      const double maxdist     =  mxGetScalar(prhs[11]);
      const double maxt        =  mxGetScalar(prhs[12]);

      int particles = mxGetNumberOfElements(prhs[0]);


      // OUTPUTS
      plhs[0]         = mxCreateDoubleMatrix(particles, 1024, mxREAL);
      double *x_trajs = mxGetPr(plhs[0]);
      plhs[1]         = mxCreateDoubleMatrix(particles, 1024, mxREAL);
      double *y_trajs = mxGetPr(plhs[1]);
      plhs[2]         = mxCreateDoubleMatrix(particles, 1024, mxREAL);
      double *tss     = mxGetPr(plhs[2]);
      plhs[3]         = mxCreateDoubleMatrix(particles,    1, mxREAL);
      double *itss    = mxGetPr(plhs[3]);
      plhs[4]         = mxCreateDoubleMatrix(particles,    1, mxREAL);
      double *datass  = mxGetPr(plhs[4]);

      mexPrintf("Calling FORTRAN code.\n");
      ray_optics_ensemble(&particles, xxs, yys, vxxs, vyys, \
                    potential_maps, voltagess,\
                    dimensions, &n_electrodes, ms, qs, &din, &maxdist, &maxt, \
                    x_trajs, y_trajs, tss, itss, datass);
      mexPrintf("FORTRAN exited successfully.\n");
}
