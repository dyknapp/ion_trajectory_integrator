#include "mex.hpp"
#include "mexAdapter.hpp"
#include "trajectory_integration_module.h"

#define EXPECTED_INPUTS 1

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        checkArguments(outputs, inputs);
        const double input = inputs[0][0];
        
        // Get array factory
        ArrayFactory factory;

        const double result = minimal_test(&input);

        TypedArray<double> output = factory.createScalar(result);
        
        outputs[0] = output;
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
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("First input must be scalar double") }));
        }

        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
};