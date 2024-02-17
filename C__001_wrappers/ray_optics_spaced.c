#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "ray_optics_spaced.h"

     //  subroutine ray_optics_spaced(position, velocity, sample_dist,
     // &                        potential_maps, voltages, dimensions,
     // &                        n_electrodes, m, q, din, maxdist, maxt,
     // &                        trajectory, death, its, datas)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      mxDouble *position        = mxGetDoubles(prhs[0 ]);
      mxDouble *velocity        = mxGetDoubles(prhs[1 ]);
      const double sample_dist  = mxGetScalar( prhs[2 ]);
      mxInt32  *is_electrode    = mxGetInt32s( prhs[3 ]);
      mxDouble *potential_maps  = mxGetDoubles(prhs[4 ]);
      mxDouble *voltages        = mxGetDoubles(prhs[5 ]);
      mxInt32  *dimensions      = mxGetInt32s( prhs[6 ]);
      const int n_electrodes    = mxGetScalar( prhs[7 ]);
      const double m            = mxGetScalar( prhs[8 ]);
      const double q            = mxGetScalar( prhs[9 ]);
      const double din          = mxGetScalar( prhs[10]);
      const double max_dist     = mxGetScalar( prhs[11]);
      const double maxt         = mxGetScalar( prhs[12]);

      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(1024, 3, mxREAL);
      double *trajectory = mxGetPr(plhs[0]);
      int death, its, datas;

      ray_optics_spaced(position, velocity, &sample_dist, is_electrode, \
                        potential_maps, voltages, dimensions, \
                        &n_electrodes, &m, &q, &din, &max_dist, &maxt, \
                        trajectory, &death, &its, &datas);


      plhs[1] = mxCreateDoubleScalar((double)death);
      plhs[2] = mxCreateDoubleScalar((double)its);
      plhs[3] = mxCreateDoubleScalar((double)datas);
}
