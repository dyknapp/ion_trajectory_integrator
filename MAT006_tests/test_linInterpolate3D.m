% Test the linInterpolate3D FORTRAN version implementation
%
% If you want to do this in unix (as you should) figure it out yourself.
%   -> take a peak at MAT005_mex_compatiblity/FORTRAN_to_mex.m
%
% 4.9.2023  dknapp,     wrote the script

%% Compile the mex file (Windows)

clear functions
!gfortran -Wall -Ofast -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
mex -r2018a C__002_test_wrappers/test_linInterpolate3D_fortran.c trajectory_integration_module.o -outdir FOR001_modules/
!del trajectory_integration.mod trajectory_integration_module.o

%% Testing

test_potential = normrnd(0, 1, [64 1]);

d = 0.1 + rand();

x = (1 + rand() * 3) * d;
y = (1 + rand() * 3) * d;
z = (1 + rand() * 3) * d;

test_linInterpolate3D_fortran(test_potential, x, y, z, d)
linInterpolate3D(reshape(test_potential, [4 4 4]), x, y, z, d)