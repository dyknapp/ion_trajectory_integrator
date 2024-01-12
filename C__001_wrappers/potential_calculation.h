void iterate_laplace (const double *guess, const double *bc_mask, const double *threshold, const double *maxits, const int *d1, const int *d2, double *its, double *refined);
void nn_interpolate (const double *array, const int *old_samples, const int *new_samples, double *new_array);
void nn_interpolate2D (const double *array, const int *odim1, const int *odim2, const int *ndim1, const int *ndim2, double *new_array);
void refined_laplace (const double *guess, const double *bc_mask, const double *threshold, const double *maxits, const int *d1, const int *d2, const int *refinements, double *refined);
