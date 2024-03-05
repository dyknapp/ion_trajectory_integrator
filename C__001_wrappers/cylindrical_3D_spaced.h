#include <stddef.h>
#ifdef __cplusplus
#include <complex>
#define __GFORTRAN_FLOAT_COMPLEX std::complex<float>
#define __GFORTRAN_DOUBLE_COMPLEX std::complex<double>
#define __GFORTRAN_LONG_DOUBLE_COMPLEX std::complex<long double>
extern "C" {
#else
#define __GFORTRAN_FLOAT_COMPLEX float _Complex
#define __GFORTRAN_DOUBLE_COMPLEX double _Complex
#define __GFORTRAN_LONG_DOUBLE_COMPLEX long double _Complex
#endif

void cylindrical_3D_spaced (const double *position, const double *velocity, const double *sample_dist, const int *is_electrode, const double *potential_maps, const double *voltages, const int *voltage_lines, const int *dimensions, const int *n_electrodes, const double *m, const double *q, const double *din, const double *maxdist, const double *maxt, double *trajectory, int *death, int *its, int *datas);

#ifdef __cplusplus
}
#endif
