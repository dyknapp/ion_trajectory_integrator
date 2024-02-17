#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "ray_optics_spaced_ensemble.h"

     //  subroutine ray_optics_spaced(particles, position, velocity, sample_dist,
     // &                        potential_maps, voltages, dimensions,
     // &                        n_electrodes, m, q, din, maxdist, maxt,
     // &                        trajectory, death, its, datas)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      const int particles       = mxGetScalar( prhs[0 ]);
      mxDouble *positions       = mxGetDoubles(prhs[1 ]);
      mxDouble *velocitys       = mxGetDoubles(prhs[2 ]);
      const double sample_dist  = mxGetScalar( prhs[3 ]);
      mxInt32  *is_electrode    = mxGetInt32s( prhs[4 ]);
      mxDouble *potential_maps  = mxGetDoubles(prhs[5 ]);
      mxDouble *voltages        = mxGetDoubles(prhs[6 ]);
      mxInt32  *dimensions      = mxGetInt32s( prhs[7 ]);
      const int n_electrodes    = mxGetScalar( prhs[8 ]);
      const double m            = mxGetScalar( prhs[9 ]);
      const double q            = mxGetScalar( prhs[10]);
      const double din          = mxGetScalar( prhs[11]);
      const double max_dist     = mxGetScalar( prhs[12]);
      const double maxt         = mxGetScalar( prhs[13]);

      // OUTPUTS
      plhs[0] = mxCreateDoubleMatrix(particles, 1024 * 3, mxREAL);
      double *trajectories = mxGetPr(plhs[0]);
      plhs[1] = mxCreateNumericMatrix(1, particles, mxINT32_CLASS, mxREAL);
      double *deaths = mxGetPr(plhs[1]);
      plhs[2] = mxCreateNumericMatrix(1, particles, mxINT32_CLASS, mxREAL);
      double *itss = mxGetPr(plhs[2]);
      plhs[3] = mxCreateNumericMatrix(1, particles, mxINT32_CLASS, mxREAL);
      double *datass = mxGetPr(plhs[3]);

      ray_optics_spaced_ensemble(
            &particles, positions, velocitys, &sample_dist, is_electrode, \
            potential_maps, voltages, dimensions, \
            &n_electrodes, &m, &q, &din, &max_dist, &maxt, \
            trajectories, deaths, itss, datass);
}
