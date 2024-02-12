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

void inf_trap_bg_gas (const double *position, const double *velocity, const double *m, const double *q, const double *coll_rate, const double *room_temp, const double *omega, const double *depth, const double *r, const double *max_t, const double *max_dist, const int *averaging_window, const int *pts, double *four_trajectory, int *its, int *rpts, int *collisions);
double legendre_f (const int *l, const int *m, const double *x);
void maxwell_distribution_test (const double *room_temp, const int *samples, double *v);
void Ylm_f (const int *l, const int *m, const double *phi, const double *theta, double *r, double *c);

#ifdef __cplusplus
}
#endif
