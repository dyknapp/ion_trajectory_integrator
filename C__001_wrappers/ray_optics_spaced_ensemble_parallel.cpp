#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
    double sm = 0;
    const TypedArray<double> inArray = inputs[0];
    for (auto& elem : inArray) {
        sm += elem;
    }
    outputs[0] = factory.createScalar(sm);
    }
};