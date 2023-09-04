//>> !gfortran -Wall -Ofast -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
//>> mex -r2018a C__002_test_wrappers/test_linInterpolate3D_fortran.c trajectory_integration_module.o -outdir FOR001_modules/
//>> !del trajectory_integration.mod trajectory_integration_module.o

#include "mex.h"
#include "matrix.h"
#include "../C__001_wrappers/trajectory_integration_module.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxDouble *potential = mxGetDoubles(prhs[0]);
	const double x = mxGetScalar(prhs[1]);
    const double y = mxGetScalar(prhs[2]);
    const double z = mxGetScalar(prhs[3]);
    const double d = mxGetScalar(prhs[4]);

	double result = lininterpolate3D(potential, &x, &y, &z, &d);

    plhs[0] = mxCreateDoubleScalar(result);
}