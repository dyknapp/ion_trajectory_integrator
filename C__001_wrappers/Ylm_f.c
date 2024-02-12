#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "Ylm_f.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      char str[1024];
      if(nrhs != 4){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "4 inputs required.");
      }
      for(int i = 0; i <= 1;  i++){
            if(!mxIsInt32(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected int32)", i);
                  mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", str);
            }
      }
      for(int i = 2; i <= 3;  i++){
            if(!mxIsDouble(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected double)", i);
                  mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", str);
            }
      }

      // INPUTS
      const int                  l = mxGetScalar(prhs[0]);
      const int                  m = mxGetScalar(prhs[1]);
      const double           theta = mxGetScalar(prhs[2]);
      const double             phi = mxGetScalar(prhs[3]);

      if(abs(m) > abs(l)){
            mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", "m must fulfill m<=l");
      }

      double r, c;
      Ylm_f(&l, &m, &phi, &theta, &r, &c);

      plhs[0] = mxCreateDoubleScalar(r);
      plhs[1] = mxCreateDoubleScalar(c);
}
