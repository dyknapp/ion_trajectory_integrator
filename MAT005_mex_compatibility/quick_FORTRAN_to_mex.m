function quick_FORTRAN_to_mex(name, setvars_path, debug, init)
    arguments
        name (1, 1) string
        setvars_path (1, 1) string = "default"
        debug (1, 1) logical = false
        init  (1, 1) logical = true
    end
    key = compilation_key();
    combination = key(name);
    FORTRAN_to_mex(combination{1}, combination{2}, "default", setvars_path, debug, init)
end