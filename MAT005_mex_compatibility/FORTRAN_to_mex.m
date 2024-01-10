% Script for compiling fortran files into .mex functions automatically
%
% dknapp,  1.09.2023: wrote script
% dknapp,  2.09.2023: added windows compatibility
% dknapp, 23.10.2023: Modifications to integrate into more sophisticated automation
%   
% To save time on MEX compilation, specify which to compile as a list in
% the parameter.

% mex -r2018a C__001_wrappers/trajectory_integration_module.c           "trajectory_integration_module"
% mex -r2018a C__001_wrappers/ensemble_trajectory_integration_module.c  "ensemble_trajectory_integration_module"
% mex -r2018a C__001_wrappers/fly_aqs.c                                 "fly_aqs"
% mex -r2018a C__001_wrappers/integrate_trajectory_lite.c               "integrate_trajectory_lite"
% mex -r2018a C__001_wrappers/iterate_laplace.c                         "iterate_laplace"
% mex -r2018a C__001_wrappers/refined_laplace.c                         "refined_laplace"

function FORTRAN_to_mex(which, setvars_path)
    arguments
        which (1,:) string = []
        setvars_path (1,1) string = "default"
    end
    % Make sure that the mex files are not in use
    clear functions

    % If which is empty, compile all MEX files
    if isempty(which)
        all = true;
        % initialize which as empty array
        which = ([" "]);
    else
        all = false;
    end
    
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
        if ~strcmp(setvars_path, "default")
            result = system(strcat(setvars_path, " intel64 vs2022"));
        else
            fprintf("*** WARNING *** Didn't detect a path to Intel's script for setting environment variables.  Trying with default... \n")
            result = system('"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022');
        end
        if result ~= 0
            fprintf("*** FAILURE *** oneAPI environment initialization \n")
        end
        !ifort /fpp /c /dll /O3 /Qparallel /Qopenmp .\FOR001_modules\utils.f                                   -o .\FOR001_modules\utils.o
        !ifort /fpp /c /dll /O3 /Qparallel /Qopenmp .\FOR001_modules\trajectory_integration_module.f           -o .\FOR001_modules\tim.o
        !ifort /fpp /c /dll /O3 /Qparallel /Qopenmp .\FOR001_modules\potential_calculation.f                   -o .\FOR001_modules\pc.o
        !ifort /fpp /c /dll /O3 /Qparallel /Qopenmp .\FOR001_modules\ion_optics.f                              -o .\FOR001_modules\io.o

        % Generate header file
        !del *.mod
        !gfortran -c .\FOR001_modules\utils.f
        !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\trajectory_integration_module.f > C__001_wrappers\trajectory_integration_module.h
        !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\potential_calculation.f         > C__001_wrappers\potential_calculation.h
        !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\ion_optics.f                    > C__001_wrappers\ion_optics.h
        !gfortran -fc-prototypes -fsyntax-only -cpp FOR001_modules\utils.f                         > C__001_wrappers\utils.h
        
        % Compile C wrapper to mex
        name = "trajectory_integration_module";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update')
            mex -r2018a C__001_wrappers/trajectory_integration_module.c          FOR001_modules/tim.o -outdir FOR001_modules
        end

        name = "ensemble_trajectory_integration_module";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update')
            mex -r2018a C__001_wrappers/ensemble_trajectory_integration_module.c FOR001_modules/tim.o -outdir FOR001_modules
        end

        name = "fly_aqs";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update')
            mex -r2018a C__001_wrappers/fly_aqs.c                                FOR001_modules/tim.o -outdir FOR001_modules
        end

        name = "integrate_trajectory_lite";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update')
            mex -r2018a C__001_wrappers/integrate_trajectory_lite.c              FOR001_modules/tim.o -outdir FOR001_modules
        end

        name = "iterate_laplace";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update')
            mex -r2018a C__001_wrappers/iterate_laplace.c                        FOR001_modules/pc.o -outdir FOR001_modules 
        end

        name = "refined_laplace";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/refined_laplace.c                        FOR001_modules/pc.o -outdir FOR001_modules
        end

        name = "ray_optics_2D";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/ray_optics_2D.c                          FOR001_modules/io.o -outdir FOR001_modules
        end

        name = "ray_optics_ensemble";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/ray_optics_ensemble.c                    FOR001_modules/io.o -outdir FOR001_modules
        end

        name = "interpolate_monotonic_1d";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/interpolate_monotonic_1d.c                FOR001_modules/utils.o -outdir FOR001_modules
        end

        name = "linspace_fortran";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/linspace_fortran.c                        FOR001_modules/utils.o -outdir FOR001_modules
        end

        name = "fly_ensemble";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/fly_ensemble.c                            FOR001_modules/tim.o -outdir FOR001_modules
        end

        name = "fly_cloud";
        if any(contains(which, name)) || all
            fprintf("\nCompiling MEX: %s\n", name);
            drawnow('update') 
            mex -r2018a C__001_wrappers/fly_cloud.c                               FOR001_modules/tim.o -outdir FOR001_modules
        end

        % Remove extra files
        !del *.o
    end
end