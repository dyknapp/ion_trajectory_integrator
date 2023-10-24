#include "mex.hpp"
#include "mexAdapter.hpp"
#include "trajectory_integration_module.h"
#include <string>

#define EXPECTED_INPUTS 3
#define EXPECTED_OUTPUTS 1

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // Get array factory
        ArrayFactory factory;

        try{
            //  subroutine time_interpolate(xs, ts, nin, points, nout,
            // &                            interpolated)
            checkArguments(outputs, inputs);

            // getNumberOfElements needs to come BEFORE;
            size_t in_size = inputs[0].getNumberOfElements();
            double *xs      = TypedArray<double>(std::move(inputs[0])).release().get();
            double *ts      = TypedArray<double>(std::move(inputs[1])).release().get();
            size_t out_size = inputs[2].getNumberOfElements();
            double *points  = TypedArray<double>(std::move(inputs[2])).release().get();

            buffer_ptr_t<double> interpolated_buffer = factory.createBuffer<double>(out_size);
            double *interpolated = interpolated_buffer.get();

            int nin  = (int) in_size;
            int nout = (int)out_size;
            time_interpolate(xs, ts, &nin, points, &nout, interpolated);

            outputs[0] = factory.createArrayFromBuffer({out_size}, std::move(interpolated_buffer));
        }
        catch(const std::runtime_error& ex){
            std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr 
                = getEngine();
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar(ex.what()) }));
        }
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr 
            = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check number of inputs
        if (inputs.size() != EXPECTED_INPUTS)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({factory.createScalar("Invalid number of inputs")}));
        }

        // Check input types and sizes
        if (   inputs[0].getType() != ArrayType::DOUBLE 
            || inputs[1].getType() != ArrayType::DOUBLE
            || inputs[2].getType() != ArrayType::DOUBLE)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Inputs must be scalar double arrays") }));
        }

        if (inputs[0].getNumberOfElements() != inputs[1].getNumberOfElements())
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input vector dimensions are inconsistent") }));
        }

        // Check number of outputs
        if (outputs.size() > EXPECTED_OUTPUTS) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Incorrect number of outputs") }));
        }
    }
};