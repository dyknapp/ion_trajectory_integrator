#define STINTLENGTH 65536
#define EXPECTED_INPUTS 12
#define EXPECTED_OUTPUTS 4

#include "nbody4.h"

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

    // Create an output stream
    std::ostringstream stream;

public:
    // The constructor
    MexFunction() {
        // mexLock();
    }

    // The destructor
    ~MexFunction() {
        // mexUnlock();
    }

    void operator()(ArgumentList outputs, ArgumentList inputs) {
        checkArguments(outputs, inputs);

        const int     particles                         = inputs[ 0][0];
        const double* positions      = getDataPtr<double>(inputs[ 1]);
        const double* velocities     = getDataPtr<double>(inputs[ 2]);
        const double* ms             = getDataPtr<double>(inputs[ 3]);
        const double* qs             = getDataPtr<double>(inputs[ 4]);
        const double  omega                             = inputs[ 5][0];
        const double  depth                             = inputs[ 6][0];
        const double  R                                 = inputs[ 7][0];
        const double  max_t                             = inputs[ 8][0];
        const double  max_dist                          = inputs[ 9][0];
        const double  record_step                       = inputs[10][0];
        const double  burst_time                        = inputs[11][0];

        double *trajectories \
            = (double*)malloc(sizeof(double)*particles*STINTLENGTH*6);
        double *times        \
            = (double*)malloc(sizeof(double)*STINTLENGTH);

        int its = 0;
        int recorded = 0;

        // double test = positions[2];
        // for(int i = 0; i < 20; i++){
        //     test = test * 10;
        //     stream << floor(test);
        //     test = test - (double)floor(test);
        // }
        // stream << "\n";
        // displayOnMATLAB(stream);

        nbody4( \
              &particles,  \
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
              &burst_time, \
              trajectories,  \
              times, \
              &its, \
              &recorded \
              );

        TypedArray<double> trajectories_out = \
            factory.createArray({(uint64_t)STINTLENGTH, \
                                 (uint64_t)particles, \
                                 (uint64_t)6}, \
                trajectories, \
                trajectories + STINTLENGTH*particles*6);
        free(trajectories);

        TypedArray<double> times_out = \
            factory.createArray({(uint64_t)STINTLENGTH}, \
                times, \
                times + STINTLENGTH);
        free(times);

        outputs[0] = trajectories_out;
        outputs[1] = times_out;
        outputs[2] = factory.createScalar(its);
        outputs[3] = factory.createScalar(recorded);
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Check number of inputs
        if (inputs.size() != EXPECTED_INPUTS)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({factory.createScalar("Invalid number of inputs")}));
        }

        // Check number of outputs
        if (outputs.size() > EXPECTED_OUTPUTS) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Invalid number of outputs") }));
        }
    }

    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};
