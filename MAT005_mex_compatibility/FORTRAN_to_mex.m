% Script for compiling fortran files into .mex functions automatically
%
% dknapp,  1.09.2023: wrote script
% dknapp,  2.09.2023: added windows compatibility
% dknapp, 23.10.2023: Modifications to integrate into more sophisticated automation
% dknapp, 13.01.2023: Reorganized to make it more general
%   
% To save time on MEX compilation, specify which to compile as a list in
% the parameter.

function FORTRAN_to_mex(fortran_file_path, ...
                        c_wrapper_path, ...
                        mex_name, ...
                        setvars_path, ...
                        debug, ...
                        load_variables)
    arguments
        fortran_file_path (1, 1) string
        c_wrapper_path    (1, 1) string
        mex_name          (1, 1) string  = "default"
        setvars_path      (1, 1) string  = "default"
        debug             (1, 1) logical = true
        load_variables    (1, 1) logical = true
    end
    % Make sure that the mex files are not in use
    clear functions
    
    % Check validity of FORTRAN file path
    if ~exist(fortran_file_path, "file")
        error("Automated compilation system: Could not find FORTRAN file.\n");
    else
        fortran_file_path = which(fortran_file_path);
    end

    % Validity of C wrapper path
    if ~exist(c_wrapper_path, "file")
        error("Automated compilation system: Could not find C wrapper.\n")
    else
        c_wrapper_path = which(c_wrapper_path);
    end
    
    if isunix
        if debug
            fprintf("Specific debug options currently unsupported for UNIX.")
        end

        [~, name, ~] = fileparts(c_wrapper_path);

        if strcmp(mex_name, "default")
            mex_name = name;
        end

        % Compile FORTRAN
        if contains(system_info.constraints, "gfortran only")
            result = system(sprintf("gfortran -Wall -Ofast -fPIC -c -cpp %s -o %s.o", fortran_file_path, name));
        else
            result = system(sprintf("source %s intel64 && ifort -fpp -fPIC -O3 -qmkl -c -parallel -qopenmp %s -o %s.o", setvars_path, fortran_file_path, name));
        end
        if result ~= 0
            error("Automated compilation system: FORTRAN compilation error.\n")
        end

        % Generate header file
        result = system(sprintf("gfortran -fc-prototypes -fPIC -fsyntax-only -cpp %s > %s.h", fortran_file_path, name));
        if result ~= 0
            error("Automated compilation system: FORTRAN compilation error.\n")
        end

        % Compile C wrapper to mex
        if contains(system_info.constraints, "gfortran only")
            eval(sprintf("mex -r2018a -lgfortran %s %s.o", c_wrapper_path, name));
        else
            eval(sprintf("mex -r2018a %s %s.o", c_wrapper_path, name));
        end

        % Delete unneeded files
        system(sprintf("rm %s.o", name));
        !rm *.mod
    elseif ispc     % For windows
        [~, name, ~] = fileparts(c_wrapper_path);

        if strcmp(mex_name, "default")
            mex_name = name;
        end

        % Generate header file
        result = system(sprintf("gfortran -fc-prototypes -fsyntax-only -cpp " + ...
                                "%s " + ...
                                "> C__001_wrappers/%s.h", fortran_file_path, name));
        if result ~= 0
            error("Automated compilation system: error in gfortran header file generation.")
        end

        % Compile to .o file
        if load_variables
            if ~strcmp(setvars_path, "default")
                result = system(strcat(setvars_path, " intel64 vs2022"));
            else
                fprintf("*** WARNING *** Didn't detect a path to Intel's script for setting environment variables.  Trying with default... \n")
                result = system('"C:/Program Files (x86)/Intel/oneAPI/setvars.bat" intel64 vs2022');
            end
            if result ~= 0
                error("Automated compilation system: oneAPI environment initialization \n")
            end
        end
        
        lines = readlines(fortran_file_path);
        if debug
            if contains(lines(1), "#define DEBUG")
                lines(1) = "#define DEBUG 1";
                writelines(lines, fortran_file_path);
            end
            result = system(sprintf("ifort /fpp /c /dll /O3 /Qparallel /check:bounds " + ...
                                    "%s -o ./FOR001_modules/%s.o", fortran_file_path, name));
        else
            if contains(lines(1), "#define DEBUG")
                lines(1) = "#define DEBUG 0";
                writelines(lines, fortran_file_path);
            end
            result = system(sprintf("ifort /fpp /c /dll /O3 /Qparallel " + ...
                                    "%s -o ./FOR001_modules/%s.o", fortran_file_path, name));
        end
        if result ~= 0
            error("Automated compilation system: error in FORTRAN file compilation.")
        end
        
        % Compile C wrapper to mex
        fprintf("\nCompiling MEX: %s\n", name);
        drawnow('update')
        mex('-r2018a', c_wrapper_path, sprintf("FOR001_modules/%s.o", name), ...
            '-outdir', 'FOR001_modules', '-output', mex_name)

        % Remove extra files
        !del /S *.o
        !del /S *.mod
    end
end
