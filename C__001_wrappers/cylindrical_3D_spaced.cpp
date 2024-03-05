#include <math.h>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "cylindrical_3D_spaced.h"

using namespace matlab::data;
using matlab::mex::ArgumentList;

 //  subroutine cylindrical_3D_spaced(position, velocity, sample_dist,
 // &            is_electrode, potential_maps, voltages, voltage_lines,
 // &            dimensions, n_electrodes, m, q, din, maxdist, maxt,
 // &            trajectory, death, its, datas)

template<typename T> const T* getDataPtr(Array arr){
    const TypedArray<T> arr_t = arr;
    matlab::data::TypedIterator<const T> it(arr_t.begin());
    return it.operator->();
}

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;

     // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Create an output stream
    std::ostringstream stream;
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {


    const double* position       = getDataPtr<double>(inputs[ 0]);
    const double* velocity       = getDataPtr<double>(inputs[ 1]);
    const double  sample_dist                       = inputs[ 2][0];
    const int*    is_electrode      = getDataPtr<int>(inputs[ 3]);
    const double* potential_maps = getDataPtr<double>(inputs[ 4]);
    const double* voltages       = getDataPtr<double>(inputs[ 5]);
    const int*    voltage_lines     = getDataPtr<int>(inputs[ 3]);
    const int*    dimensions        = getDataPtr<int>(inputs[ 7]);
    const int     n_electrodes                      = inputs[ 8][0];
    const double  m                                 = inputs[ 9][0];
    const double  q                                 = inputs[10][0];
    const double  din                               = inputs[11][0];
    const double  maxdist                           = inputs[12][0];
    const double  maxt                              = inputs[13][0];

    // size(trajectories) = [1000 1024 3]
    double *trajectory = (double*)malloc(sizeof(double) * 1024 * 4);
    int    death, its, datas;

    cylindrical_3D_spaced(position, velocity, &sample_dist, \
                          is_electrode, potential_maps, voltages, voltage_lines, \
                          dimensions, &n_electrodes, &m, &q, &din, &maxdist, &maxt, \
                          trajectory, &death, &its, &datas);

    TypedArray<double> trajectory_out = \
        factory.createArray({(uint64_t)1024, (uint64_t)4}, \
            trajectory, \
            trajectory + 1024 * 4);

    outputs[0] = trajectory_out;
    outputs[1] = factory.createScalar<int>(death);
    outputs[2] = factory.createScalar<int>(its);
    outputs[3] = factory.createScalar<int>(datas);

    free(trajectory);
    }

    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};