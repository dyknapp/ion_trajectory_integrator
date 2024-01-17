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
            "Not all environment checks were passed.  Compilation will be attempted, but is unlikely to succeed.")
else
    setvars_path = "default";
end

key = compilation_key();
keys = key.keys;
for n = 1:length(keys)
    combination = key(keys{n});
    FORTRAN_to_mex(combination{1}, combination{2}, ...
                   "default", setvars_path, false, (n == 1))
end

fprintf("\n\n\nREMINDERS for before you get to work:\n" + ...
        "*) Remember to CLOSE MATLAB before pushing to github :-)\n" + ...
        "  (.gitignore is automatically updated on close.)\n" + ...
        "\n" + ...
        "*) All MATLAB code in this repository is designed to only load\n" + ...
        "SIMION potentials if they are not already loaded.  If you make changes\n" + ...
        "to the potential itself or the code that loads it, make sure to\n" + ...
        "clear variables so that the potential is reloaded from scratch :)\n");
