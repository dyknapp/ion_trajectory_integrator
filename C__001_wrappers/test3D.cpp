#include "mex.hpp"
#include "mexAdapter.hpp"
#include "trajectory_integration_module.h"
#include <string>
#include <math.h>

#define EXPECTED_INPUTS  1
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

            double *a                  = TypedArray<double>(std::move(inputs[0 ])).release().get();

            size_t out_size = 10000;

            buffer_ptr_t<double> b_buffer  = factory.createBuffer<double>(out_size);

            double *b_pass = b_buffer.get();

            test3D(a, b_pass);

            outputs[0] = factory.createArrayFromBuffer({out_size}, std::move(b_buffer));
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

        // Check number of outputs
        if (outputs.size() > EXPECTED_OUTPUTS) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Incorrect number of outputs") }));
        }
    }

    void displayOnMATLAB(std::ostringstream &stream){
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr 
            = getEngine();
        ArrayFactory factory;
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({factory.createScalar(stream.str())}));
        stream.str("");
    }
};