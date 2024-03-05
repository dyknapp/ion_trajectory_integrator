function [ fix ] = langevinBath( temperature, dampingTime, group, seed )
%LANGEVINBATH Creates a langevin bath of a given temperature (in Kelvin).
% The langevin bath applies a viscous damping force to each atom
% proportional to its velocity plus a stochastic, white noise force. The
% magnitude of this noise is such that at times significantly longer than
% dampingTime the system will thermalise to the specified temperature. The
% damping time (s) is the time taken for velocity to relax to 1/e its initial
% value in a zero temperature bath.
%
% group is an optional parameter that allows the Langevin bath to be
% coupled to only some atoms. If left unspecified the bath is coupled to
% all atoms.
%
% Syntax: langevinBath( temperature, dampingTime, group )
%
% Example: sim.Add( langevinBath(4, 1e-4, sim.Group(calciumIons)) )
%
% See Also: laserCool, http://lammps.sandia.gov/doc/fix_langevin.html

fix = LAMMPSFix();
fix.createInputFileText = @langevinBathCfg;

if nargin < 3
    fix.InputFileArgs = { temperature, dampingTime };
elseif nargin < 4    
    fix.InputFileArgs = { temperature, dampingTime, group };
else
    fix.InputFileArgs = { temperature, dampingTime, group, seed };
end


end

