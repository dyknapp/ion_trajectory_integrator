% Script for compiling fortran files into .mex functions automatically
%
% dknapp, 1.9.2023: wrote script
% dknapp, 2.9.2023: added windows compatibility

% Make sure that the mex files are not in use
clear functions

if isunix       % For unix-based: NEEDS UPDATING
    % % Compile FORTRAN
    % !sudo gfortran -Wall -Ofast -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
    % % Generate header file
    % !sudo gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h
    % % Compile C wrapper to mex
    % mex -r2018a -lgfortran C__001_wrappers/trajectory_integration_module.c trajectory_integration_module.o -outdir FOR001_modules/
    % % Delete unneeded files
    % !rm trajectory_integration.mod trajectory_integration_module.o
elseif ispc     % For windows
    % Generate header file
    !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\trajectory_integration_module.f > C__001_wrappers\trajectory_integration_module.h
    !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\potential_calculation.f         > C__001_wrappers\potential_calculation.h
    % Compile to .o file
    
    !"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022
    !ifort /fpp /c /dll /O3 /Qopenmp /Qparallel .\FOR001_modules\trajectory_integration_module.f -o .\FOR001_modules\tim.o
    !ifort /fpp /c /dll /O3 /Qopenmp /Qparallel .\FOR001_modules\potential_calculation.f         -o .\FOR001_modules\pc.o
    % Compile C wrapper to mex
    mex -r2018a C__001_wrappers/trajectory_integration_module.c          FOR001_modules/tim.o -outdir FOR001_modules
    % mex -r2018a C__001_wrappers/ensemble_trajectory_integration_module.c FOR001_modules/tim.o -outdir FOR001_modules
    % mex -r2018a C__001_wrappers/iterate_laplace.c                        FOR001_modules/pc.o -outdir FOR001_modules
    % mex -r2018a C__001_wrappers/nn_interpolate.c                         FOR001_modules/pc.o -outdir FOR001_modules
    % mex -r2018a C__001_wrappers/nn_interpolate2D.c                       FOR001_modules/pc.o -outdir FOR001_modules
    % mex -r2018a C__001_wrappers/refined_laplace.c                        FOR001_modules/pc.o -outdir FOR001_modules
    % Remove extra files
    !del FOR001_modules\tim.o
    !del FOR001_modules\pc.o
    !del trajectory_integration.mod
    !del potential_calculation.mod
end
