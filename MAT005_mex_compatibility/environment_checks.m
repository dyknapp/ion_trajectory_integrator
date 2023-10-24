function [all_checks_passed, setvars_path] = environment_checks()
      all_checks_passed = true;
    
    if ispc
        fprintf("Detected PC environment.  Automatic system is compatible.\n")
    else
        all_checks_passed = false;
    end
    
    
    % Check C compiler for MEX
    fprintf("C compiler:                ")
    myCCompiler = mex.getCompilerConfigurations('C','Selected');
    if     strcmp(myCCompiler.Manufacturer, 'Intel') ...
        && str2double(myCCompiler.Version) >= 23.0
        fprintf("PASSED\n")
    else
        fprintf("*** FAILED ***\n")
        all_checks_passed = false;
        fprintf("See the README.\n" + ...
                "For a MEX file to be compiled, MATLAB needs a compiler specified.\n" + ...
                "This project is designed to run using:\n" + ...
                "   Intel oneAPI 2023 with Microsoft Visual Studio 2022 (C)\n" + ...
                "Please install compiler and/or configure MATLAB.\n" + ...
                "Otherwise, you need your own solution for compiling the MEX files.\n")
    end
    
    % Check path for C compiler environment variables
    setvars_path = strcat(extractBefore(myCCompiler.Location, '\\'), "\setvars.bat");
    fprintf("Environment variables:     ")
    if     exist(setvars_path, 'file') > 0
        fprintf("PASSED\n")
        setvars_path = strcat('"', setvars_path, '"');
    else
        fprintf("*** FAILED ***\n")
        all_checks_passed = false;
        fprintf("Your Intel compiler installation might be broken.\n")
    end

    % Check gfortran
    fprintf("gFortran:                  ")
    [result1, output1] = system("gfortran --version");
    if     result1 == 0
        % fprintf("PASSED\n%s\n", output1)
        fprintf("PASSED\n")
    else
        fprintf("*** FAILED ***\n")
        all_checks_passed = false;
        fprintf("See the README\n" + ...
                "Although Intel OneAPI FORTRAN compiler is used for compilation,\n" + ...
                "use of the automatic MEX compilation in MATLAB for this project\n" + ...
                "requires gfortran for generating header files.  Install and add to PATH.")
    end

    % Check ifort
    fprintf("Intel FORTRAN compiler:    ")
    [result2, output2] = system("ifort /QV");
    if     result2 == 0
        % fprintf("PASSED\n%s\n", output2)
        fprintf("PASSED\n")
    else
        fprintf("*** FAILED ***\n")
        all_checks_passed = false;
        fprintf("See the README\n" + ...
                "Intel OneAPI FORTRAN compiler is used for MEX compilation.  Install it.")
    end

    fprintf('\n');
end