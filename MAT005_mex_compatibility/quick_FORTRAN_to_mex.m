function quick_FORTRAN_to_mex(name)
    key = compilation_key();
    combination = key(name);
    FORTRAN_to_mex(combination{1}, combination{2})
end