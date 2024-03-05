#include <math.h>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "ion_optics.h"

using namespace matlab::data;
using matlab::mex::ArgumentList;

// subroutine ray_optics_spaced_ensmble(particles, 
// &            positions, velocities,
// &            sample_dist, is_electrode, potential_maps, voltages, 
// &            dimensions, n_electrodes, m, q, din, maxdist, maxt,
// &            trajectories, deaths, itss, datass)

template<typename T> const T* getDataPtr(Array arr){
    const TypedArray<T> arr_t = arr;
    matlab::data::TypedIterator<const T> it(arr_t.begin());
    return it.operator->();
}

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {

    const int     particles                         = inputs[ 0][0];
    const double* positions      = getDataPtr<double>(inputs[ 1]);
    const double* velocities     = getDataPtr<double>(inputs[ 2]);
    const double  sample_dist                       = inputs[ 3][0];
    const int*    is_electrode      = getDataPtr<int>(inputs[ 4]);
    const double* potential_maps = getDataPtr<double>(inputs[ 5]);
    const double* voltages       = getDataPtr<double>(inputs[ 6]);
    const int*    dimensions        = getDataPtr<int>(inputs[ 7]);
    const int     n_electrodes                      = inputs[ 8][0];
    const double  m                                 = inputs[ 9][0];
    const double  q                                 = inputs[10][0];
    const double  din                               = inputs[11][0];
    const double  maxdist                           = inputs[12][0];
    const double  maxt                              = inputs[13][0];

    // size(trajectories) = [1000 1024 3]
    double *trajectories = (double*)malloc(sizeof(double) * particles * 1024 * 3);
    int    *deaths       = (   int*)malloc(sizeof(int) * particles);
    int    *itss         = (   int*)malloc(sizeof(int) * particles);
    int    *datass       = (   int*)malloc(sizeof(int) * particles);

    // Temp variables for loop
    #pragma omp parallel for \
        shared(positions, velocities, sample_dist, is_electrode, potential_maps, \
               voltages, dimensions, n_electrodes, m, q, din, maxdist, maxt, \
               trajectories, deaths, itss, datass, particles) \
        default(none)
    for(int idx = 0; idx < particles; idx++){
        
        double* position_temp   = new double[2];
        double* velocity_temp   = new double[2];
        double* trajectory_temp = new double[3 * 1024];
        int death, its, datas;

        // // Initialize for current iteration
        for(int coord = 0; coord < 2; coord++){
            position_temp[coord] =  positions[idx*2 + coord];
            velocity_temp[coord] = velocities[idx*2 + coord];
        }

        ray_optics_spaced( \
            position_temp, velocity_temp, \
            &sample_dist, is_electrode, \
            potential_maps, voltages, dimensions, \
            &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
            trajectory_temp, &death, &its, &datas);

        // Write results onto main arrays
        for(int jdx = 0; jdx < 1024; jdx++){
            for(int coord = 0; coord < 3; coord++){
                trajectories[idx + jdx*1000 + coord*1000*1024] \
                    = trajectory_temp[jdx + coord*1024];
                trajectory_temp[jdx + coord*1024] = 0;
            }
        }
        deaths[idx] = death;
        itss[idx] = its;
        datass[idx] = datas;
    }


    TypedArray<double> trajectories_out = \
        factory.createArray({(uint64_t)particles, (uint64_t)1024, (uint64_t)3}, \
            trajectories, \
            trajectories + particles * 1024 * 3);
    free(trajectories);

    TypedArray<int> deaths_out = \
        factory.createArray({(uint64_t)1, (uint64_t)particles},\
            deaths,\
            deaths + particles);
    free(deaths);

    TypedArray<int> itss_out = \
        factory.createArray({(uint64_t)1, (uint64_t)particles},\
            itss,\
            itss + particles);
    free(itss);

    TypedArray<int> datass_out = \
        factory.createArray({(uint64_t)1, (uint64_t)particles},\
            datass,\
            datass + particles);
    free(datass);

    outputs[0] = trajectories_out;
    outputs[1] = deaths_out;
    outputs[2] = itss_out;
    outputs[3] = datass_out;
    }
};