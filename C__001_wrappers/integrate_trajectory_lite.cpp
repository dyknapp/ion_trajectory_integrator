#include "mex.hpp"
#include "mexAdapter.hpp"
#include "trajectory_integration_module.h"
#include <string>
#include <math.h>

#define EXPECTED_INPUTS 19
#define EXPECTED_OUTPUTS 6

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

            double *xx                 = TypedArray<double>(std::move(inputs[0 ])).release().get();
            double *yy                 = TypedArray<double>(std::move(inputs[1 ])).release().get();
            double *zz                 = TypedArray<double>(std::move(inputs[2 ])).release().get();
            double *vxx                = TypedArray<double>(std::move(inputs[3 ])).release().get();
            double *vyy                = TypedArray<double>(std::move(inputs[4 ])).release().get();
            double *vzz                = TypedArray<double>(std::move(inputs[5 ])).release().get();
            double *potential_maps     = TypedArray<double>(std::move(inputs[6 ])).release().get();
            double *voltages           = TypedArray<double>(std::move(inputs[7 ])).release().get();
            double *step_times_in      = TypedArray<double>(std::move(inputs[8 ])).release().get();
            int     time_steps         =                              inputs[9 ][0];
            int    *dimensions         = TypedArray<int>   (std::move(inputs[10])).release().get();
            int    *is_electrode       = TypedArray<int>   (std::move(inputs[11])).release().get();
            int     n_electrodes       =                              inputs[12][0];
            double  m                  =                              inputs[13][0];
            double  q                  =                              inputs[14][0];
            double  din                =                              inputs[15][0];
            double  maxdist            =                              inputs[16][0];
            double  maxt               =                              inputs[17][0];
            double  data_t_interval    =                              inputs[18][0];

            size_t out_size = ceil(maxt / data_t_interval);
            int out_length = (int)out_size;

            buffer_ptr_t<double> xs_buffer  = factory.createBuffer<double>(out_size);
            buffer_ptr_t<double> ys_buffer  = factory.createBuffer<double>(out_size);
            buffer_ptr_t<double> zs_buffer  = factory.createBuffer<double>(out_size);
            buffer_ptr_t<double> ts_buffer  = factory.createBuffer<double>(out_size);
            
            double *x_traj = xs_buffer.get();
            double *y_traj = ys_buffer.get();
            double *z_traj = zs_buffer.get();
            double *ts     = ts_buffer.get();

            double its = -1.0d;
            double recorded = -1.0d;

            integrate_trajectory_lite(xx, yy, zz, vxx, vyy, vzz, 
                potential_maps, voltages, step_times_in, 
                &time_steps, dimensions, is_electrode, &n_electrodes, 
                &m, &q, &din, &maxdist, &maxt, &data_t_interval, &out_length, 
                x_traj, y_traj, z_traj, ts, &its, &recorded);

            outputs[0] = factory.createArrayFromBuffer({out_size}, std::move(xs_buffer));
            outputs[1] = factory.createArrayFromBuffer({out_size}, std::move(ys_buffer));
            outputs[2] = factory.createArrayFromBuffer({out_size}, std::move(zs_buffer));
            outputs[3] = factory.createArrayFromBuffer({out_size}, std::move(ts_buffer));

            outputs[4] = factory.createScalar(its);
            outputs[5] = factory.createScalar(recorded);
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

        // // Check input types and sizes
        // if (   inputs[0].getType() != ArrayType::DOUBLE 
        //     || inputs[1].getType() != ArrayType::DOUBLE
        //     || inputs[2].getType() != ArrayType::DOUBLE)
        // {
        //     matlabPtr->feval(u"error",
        //         0,
        //         std::vector<Array>({ factory.createScalar("Inputs must be scalar double arrays") }));
        // }

        // if (inputs[0].getNumberOfElements() != inputs[1].getNumberOfElements())
        // {
        //     matlabPtr->feval(u"error",
        //         0,
        //         std::vector<Array>({ factory.createScalar("Input vector dimensions are inconsistent") }));
        // }

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