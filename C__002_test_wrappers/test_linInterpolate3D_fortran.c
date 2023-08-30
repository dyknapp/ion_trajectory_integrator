#include <stdio.h>
#include "linInterpolate3D_fortran.h"

int main(int argc, char *argv[])
{
	double test_input[4][4][4];
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int k = 0; k < 4; k++)
			{
				// +3 to offset the indexing of MATLAB and FORTRAN, which both start at 1.
				test_input[i][j][k] = (double)((i + j + k + 3) * (i + j + k + 3));
			}
		}
	}

	const double x = 2.1;
	const double y = 1.1;
	const double z = 0.1;
	const double d = 1.0;

	double result = lininterpolate3D((double *)&test_input, &x, &y, &z, &d);
	printf("%.3f\n", result);

	return 0;
}