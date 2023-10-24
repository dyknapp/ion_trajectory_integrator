clear; clc;
addpath(genpath('.'));
rmpath('./.git')

[all_checks_passed, setvars_path] = environment_checks();
if ~all_checks_passed
    % Update the default path in FORTRAN_to_mex.m
    setvars_path = "default";
end

if ~all_checks_passed
    fprintf("*** ERROR ***\n" + ...
            "Not all environment checks were passed.  Automatic initialization procedure aborted.")
else
    FORTRAN_to_mex([], setvars_path)
end

fprintf("\n\n\nREMINDERS for before you get to work:\n" + ...
        "Remember to CLOSE MATLAB before pushing to github :-)\n" + ...
        "  (.gitignore is automatically updated on close.)\n");
