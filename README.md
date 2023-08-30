# MATLAB Ion Trajectory Integrator

This project takes potentials exported from SIMION and computes ion trajectories in MATLAB.  This enables more flexibility in running and analyzing simulations.

The core of the code is written in FORTRAN, with wrapper functions written in C.

## Operability Status
- ./C__002_test_wrappers is currently out of commission (coming soon)

## Parking space for useful commands (Assuming you are in the ion_trajectory_integrator directory):
If you modify the FORTRAN, you need to first compile it:
```
$ sudo gfortran -Wall -O3 -fPIC -c -cpp ./FOR001_modules/trajectory_integration_module.f
```
Next you need to generate a header file for the C wrapper:
```
$ sudo gfortran -fc-prototypes -fsyntax-only -cpp ./trajectory_integration_module.f > ../C__001_wrappers/trajectory_integration_module.h
```
Finally, a MATLAB command for compiling the C wrapper into a `.mex` file:
```
>> addpath(genpath('.'));
>> mex -r2018a -lgfortran C__001_wrappers/trajectory_integration_module.c FOR001_modules/trajectory_integration_module.o
```
Now, you should have a MATLAB mex function with the same name as the C file you just compiled.