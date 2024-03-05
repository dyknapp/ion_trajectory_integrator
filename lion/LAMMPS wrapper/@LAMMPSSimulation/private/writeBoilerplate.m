function writeBoilerplate( fHandle )
%WRITEBOILERPLATE Adds boiler plate text for the LIon wrapper to the config
%file
fprintf(fHandle, '# LAMMPS Input file generated using LIon, a Matlab wrapper for ion trap simulations\n');
fprintf(fHandle, '# Lion Version: %s\n', lionVersion);
fprintf(fHandle, '###################################################################################\n\n');

end

