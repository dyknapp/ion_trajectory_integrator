% Script for compiling fortran files into .mex functions automatically
%
% dknapp, 1.9.2023: wrote script
% dknapp, 2.9.2023: added windows compatibility
if isunix
    !sudo gfortran -Wall -O3 -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f
    !sudo gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h
    mex -r2018a -lgfortran C__001_wrappers/trajectory_integration_module.c FOR001_modules/trajectory_integration_module.o -outdir FOR001_modules/
    !rm trajectory_integration.mod trajectory_integration_module.o
elseif ispc
    !wsl -u root sudo gfortran -Wall -O3 -fPIC -c -cpp FOR001_modules/trajectory_integration_module.f; exit
    !wsl -u root sudo gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules/trajectory_integration_module.f > C__001_wrappers/trajectory_integration_module.h; exit
    mex -r2018a C__001_wrappers/trajectory_integration_module.c trajectory_integration_module.o -outdir FOR001_modules/
    !wsl rm trajectory_integration.mod trajectory_integration_module.o; exit
end
