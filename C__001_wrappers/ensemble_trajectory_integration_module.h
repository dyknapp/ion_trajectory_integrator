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

void field_at (const double *t, const double *voltages, const double *step_times, const int *n_electrodes, const int *time_steps, const double *x, const double *y, const double *z, const double *d, const double *potential_maps, const int *dimensions, double *ex, double *ey, double *ez);
void integrate_trajectory (const double *xx, const double *yy, const double *zz, const double *vxx, const double *vyy, const double *vzz, const double *potential_maps, const double *voltages, const double *step_times_in, const int *time_steps, const int *dimensions, const int *is_electrode, const int *n_electrodes, const double *m, const double *q, const double *din, const double *maxdist, const double *maxt, double *x_traj, double *y_traj, double *z_traj, double *ts, double *exs, double *eys, double *ezs, double *its);
void interpolate_voltages (const double *t, const double *voltages, const double *step_times, const int *n_electrodes, const int *time_steps, double *interpolated_voltages);
double lininterpolate3D (const double *matrix, const double *xin, const double *yin, const double *zin, const double *d_grid);
void potEfunc (const double *potential, const double *x, const double *y, const double *z, const double *d, double *ex, double *ey, double *ez);

#ifdef __cplusplus
}
#endif
