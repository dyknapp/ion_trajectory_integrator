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

void ray_optics_2D (const double *xx, const double *yy, const double *vxx, const double *vyy, const int *is_electrode, const double *potential_maps, const double *voltages, const int *dimensions, const int *n_electrodes, const double *m, const double *q, const double *din, const double *maxdist, const double *maxt, double *x_traj, double *y_traj, double *ts, double *its);
void ray_optics_ensemble (const int *particles, const double *xxs, const double *yys, const double *vxxs, const double *vyys, const double *potential_maps, const double *voltagess, const int *dimensions, const int *n_electrodes, const double *ms, const double *qs, const double *din, const double *maxdist, const double *maxt, double *x_trajs, double *y_trajs, double *tss, double *itss, double *datass);
void ray_optics_spaced (const double *position, const double *velocity, const double *sample_dist, const int *is_electrode, const double *potential_maps, const double *voltages, const int *dimensions, const int *n_electrodes, const double *m, const double *q, const double *din, const double *maxdist, const double *maxt, double *trajectory, int *death, int *its, int *datas);
void ray_optics_spaced_ensemble (const int *particles, const double *positions, const double *velocities, const double *sample_dist, const int *is_electrode, const double *potential_maps, const double *voltages, const int *dimensions, const int *n_electrodes, const double *m, const double *q, const double *din, const double *maxdist, const double *maxt, double *trajectories, int *deaths, int *itss, int *datass);

#ifdef __cplusplus
}
#endif
