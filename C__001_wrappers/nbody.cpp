#define STINTLENGTH 65536

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <unordered_map>
#include <fstream> 

using matlab::mex::ArgumentList;
using namespace matlab::data;

template<typename T> const T* getDataPtr(Array arr){
    const TypedArray<T> arr_t = arr;
    matlab::data::TypedIterator<const T> it(arr_t.begin());
    return it.operator->();
}

class MexFunction : public matlab::mex::Function {

    // Pointer to MATLAB engine
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    ArrayFactory factory;

public:
    // The constructor
    MexFunction() {
        mexLock();
    }

    // The destructor
    ~MexFunction() {
    }

    void operator()(ArgumentList outputs, ArgumentList inputs) {
        const int     particles                         = inputs[ 0][0];
        const double* positions      = getDataPtr<double>(inputs[ 1]);
        const double* velocities     = getDataPtr<double>(inputs[ 2]);
        const double* ms             = getDataPtr<double>(inputs[ 3]);
        const double* qs             = getDataPtr<double>(inputs[ 4]);
        const double  omega                             = inputs[ 3][0];
        const double  depth                             = inputs[ 3][0];
        const double  R                                 = inputs[ 3][0];
        const double  max_t                             = inputs[ 3][0];
        const double  max_dist                          = inputs[ 3][0];
        const double  record_step                       = inputs[ 3][0];
        const double  burst_time                        = inputs[ 3][0];

        double *trajectories = (double*)malloc(sizeof(double)*particles*STINTLENGTH*4);

        nbody( \
              particles,  \
              positions, \
              velocities,  \
              ms,  \
              qs, \
              &omega,  \
              &depth,  \
              &R,  \
              &max_t, \
              &max_dist, \
              &record_step, \
              four_trajectory,  \
              &its \
              );

        TypedArray<double> trajectories_out = \
            factory.createArray({(uint64_t)STINTLENGTH, (uint64_t)particles, (uint64_t)4}, \
                trajectories, \
                trajectories + (sizeof(trajectories) / sizeof(double)));
        free(trajectories);

        outputs[0] = trajectories_out;
        outputs[1] = factory.createScalar(its);
    }

    void displayOnMATLAB(const std::ostringstream& stream){
        matlabPtr->feval(u"fprintf", 0, 
            std::vector<Array>({ factory.createScalar(stream.str()) }));
    }
};