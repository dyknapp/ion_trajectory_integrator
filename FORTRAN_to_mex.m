% Script for compiling fortran files into .mex functions automatically
% UNIX: Assumes that gfortran and gcc are installed.
%
% dknapp, 1.9.2023: wrote script
if isunix
    !sudo gfortran -Wall -O3 -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
    !sudo gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h
    mex -r2018a -lgfortran C__001_wrappers/trajectory_integration_module.c FOR001_modules/trajectory_integration_module.o -outdir FOR001_modules/
else
    disp("Platform not supported.")
end
