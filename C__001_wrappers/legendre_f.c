#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "legendre_f.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      char str[1024];
      if(nrhs != 3){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("trajectory_integration_module:InvalidInput", "3 inputs required.");
      }
      for(int i = 0; i <= 1;  i++){
            if(!mxIsInt32(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected int32)", i);
                  mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", str);
            }
      }
      if(!mxIsDouble(prhs[2])){
            sprintf(str, "Input %d: invalid type. (expected double)", 2);
            mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", str);
      }

      // INPUTS
      const int                  l = mxGetScalar(prhs[0]);
      const int                  m = mxGetScalar(prhs[1]);
      const double               x = mxGetScalar(prhs[2]);

      if(m > l){
            mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", "m must fulfill m<=l");
      }
      if(x < -1 || x > 1){
            mexErrMsgIdAndTxt("MEX:legendre:InvalidInput", "m must fulfill -1<=x<=1");
      }

      plhs[0] = mxCreateDoubleScalar(legendre_f(&l, &m, &x));
}
