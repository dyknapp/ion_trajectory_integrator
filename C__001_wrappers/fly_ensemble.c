#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "trajectory_integration_module.h"

#define MAX_TRAJECTORY_POINTS 1048576
//      ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE FORTRAN FILE ALSO

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      char str[1024];

      if(nrhs != 20){
            mexPrintf("Received %d inputs.\n", nrhs);
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", "20 inputs required.");
      }

      if(nlhs != 5){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidOutput", "Anticipated 5 outputs.");
      }

      // INPUTS
      for(int i = 0; i <= 1;  i++){
            if(!mxIsInt32(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected int32)", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      const int    pts = mxGetScalar(prhs[0]);
      const int    n_particles = mxGetScalar(prhs[1]);
      if(n_particles > 16384){
            mexPrintf("WARNING: Large number of particles may cause instability.");
      }

      for(int i = 2; i <= 7;  i++){
            if(!mxIsDouble(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected double)", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
            if(mxGetNumberOfElements(prhs[i]) != n_particles){
                  sprintf(str, "Input %d: inconsistent length.", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      mxDouble *xx  = mxGetDoubles(prhs[2]);
      mxDouble *yy  = mxGetDoubles(prhs[3]);
      mxDouble *zz  = mxGetDoubles(prhs[4]);
      mxDouble *vxx = mxGetDoubles(prhs[5]);
      mxDouble *vyy = mxGetDoubles(prhs[6]);
      mxDouble *vzz = mxGetDoubles(prhs[7]); // 8

      for(int i = 8; i <= 10;  i++){
            if(!mxIsDouble(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected double)", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      mwSize *pm_dimensions = mxGetDimensions(prhs[8]);
      mwSize *v_dimensions = mxGetDimensions(prhs[9]);
      mwSize *sti_dimensions = mxGetDimensions(prhs[10]);

      mxDouble *potential_maps = mxGetDoubles(prhs[8]);
      mxDouble *voltages = mxGetDoubles(prhs[9]);
      mxDouble *step_times_in = mxGetDoubles(prhs[10]); //3


      for(int i = 11; i <= 14;  i++){
            if(!mxIsInt32(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected int32)", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      for(int i = 15; i <= 19;  i++){
            if(!mxIsDouble(prhs[i])){
                  sprintf(str, "Input %d: invalid type. (expected double)", i);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }

      mwSize *ie_dimensions = mxGetDimensions(prhs[13]);
      const int time_steps = mxGetScalar(prhs[11]);
      mxInt32 *dimensions = mxGetInt32s(prhs[12]);
      mxInt32 *is_electrode = mxGetInt32s(prhs[13]); // 3

      const int n_electrodes = mxGetScalar(prhs[14]);
      const double m = mxGetScalar(prhs[15]);
      const double q = mxGetScalar(prhs[16]);
      const double din = mxGetScalar(prhs[17]);
      const double maxdist = mxGetScalar(prhs[18]);
      const double maxt = mxGetScalar(prhs[19]);  // 6 (sum: 18)

      if(pm_dimensions[0] != n_electrodes){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", \
                  "Potential Maps dimension 1 is inconsistent with inputted number of electrodes.");
      }
      for(int i = 0; i < 3; i++){
            if(pm_dimensions[i + 1] != dimensions[i]){
                  sprintf(str, "potential_maps dimension %d is inconsistent with inputted dimensions.", i + 2);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      for(int i = 0; i < 3; i++){
            if(ie_dimensions[i] != dimensions[i]){
                  sprintf(str, "is_electrode dimension %d is inconsistent with inputted dimensions.", i + 2);
                  mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", str);
            }
      }
      if(v_dimensions[0] != time_steps){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", \
                  "Voltages array dimension 1 is inconsistent with inputted number of voltage timesteps.");
      }
      if(v_dimensions[1] != n_electrodes){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", \
                  "Voltages array dimension 2 is inconsistent with inputted number of electrodes.");
      }
      if(sti_dimensions[1] != time_steps){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", \
                  "Step times array (row vector) length is inconsistent with inputted number of voltage timesteps.");
      }
      if(din < 0){
            mexErrMsgIdAndTxt("MEX:fly_ensemble:InvalidInput", "d must be greater than 0.");
      }


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

      // mexPrintf("Calling FORTRAN code.\n");
      fly_ensemble(&pts, &n_particles, xx, yy, zz, vxx, vyy, vzz, \
                   potential_maps, voltages, step_times_in, \
                   &time_steps, dimensions, is_electrode, \
                   &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                   x_trajs, y_trajs, z_trajs, tss, itss);
      // mexPrintf("FORTRAN exited successfully.\n");
}
