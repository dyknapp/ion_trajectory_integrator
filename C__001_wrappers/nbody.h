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

void nbody (const int *particles, const double *positions, const double *velocities, const double *ms_in, const double *qs_in, const double *omega, const double *depth, const double *r, const double *max_t, const double *max_dist, const double *record_step, const double *burst_time, double *trajectory, double *times, int *its);

#ifdef __cplusplus
}
#endif
