function [ fix ] = fixNVEIntegrator(groupName, preamble)
%FIXNVEINTEGRATOR This fix updates the velocities and positions of atoms
%according to the forces applied by other fixes each timestep.
% Note: This fix is added automatically by the simulation and need not be
% added by the user.
% 
% Syntax: nve = fixNVEIntegrator()
%
% See Also: http://lammps.sandia.gov/doc/fix_nve.html

if nargin < 1
   groupName = 'all';
end

if nargin < 2
    preamble = '';
end

fix = LAMMPSFix();
fix.createInputFileText =  @nveIntegratorCfg;
fix.InputFileArgs = { groupName, preamble };

end

