% Script for compiling fortran files into .mex functions automatically
%
% dknapp, 1.9.2023: wrote script
% dknapp, 2.9.2023: added windows compatibility

% Make sure that the mex files are not in use
clear functions

if isunix       % For unix-based
    % Compile FORTRAN
    !sudo gfortran -Wall -Ofast -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
    % Generate header file
    !sudo gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h
    % Compile C wrapper to mex
    mex -r2018a -lgfortran C__001_wrappers/trajectory_integration_module.c trajectory_integration_module.o -outdir FOR001_modules/
    % Delete unneeded files
    !rm trajectory_integration.mod trajectory_integration_module.o
elseif ispc     % For windows
    % Compile FORTRAN
    !gfortran -Wall -Ofast -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
    % Generate header file
    !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h
    % Compile C wrapper to mex
    mex -r2018a C__001_wrappers/trajectory_integration_module.c trajectory_integration_module.o -outdir FOR001_modules/
    % Delete unneeded files
    !del trajectory_integration.mod trajectory_integration_module.o
end
