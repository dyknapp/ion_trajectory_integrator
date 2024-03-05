function [ fix ] = efield( Ex, Ey, Ez )
%EFIELD Adds a uniform, time-independent electric field to the simulation.
% This is used to implement a bias field. Ex, Ey, Ez are the
% magnitudes of the electric field specified in the cartesian basis (units
% of V/m).
% 
% SYNTAX: efield( Ex, Ey, Ez )
%
% Example: efield( 10, 0, 0 )
%
% See Also: http://lammps.sandia.gov/doc/fix_efield.html

fix = LAMMPSFix();
fix.createInputFileText = @efieldCfg;
fix.InputFileArgs = { Ex, Ey, Ez };

end