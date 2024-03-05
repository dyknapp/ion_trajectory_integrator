function [ runCmd ] = thermalVelocities( temperature, zeroTotalMom, seed)
%THERMALVELOCITIES Set the velocities of atoms to a thermal
%distribution of the specified temperature (Kelvin).
% Set zeroTotalMom to 'yes' if the resulting ensemble should have zero
% total momentum/center of mass motion.
%
% Syntax: 
%  thermalVelocities(temperature, zeroTotalMom, seed)
%  thermalVelocities(temperature, zeroTotalMom)
%  thermalVelocities(temperature)
% 
% Example:
%  thermalVelocities(300) % 300 Kelvin;
% 
% See Also: http://lammps.sandia.gov/doc/velocity.html

if nargin < 3
    seed = 1+floor(rand(1) * 10000);
end

if nargin < 2
    zeroTotalMom = 'yes';
end

if ~strcmp(zeroTotalMom,'yes') && ~strcmp(zeroTotalMom,'no')
    error('zeroTotalMom must be "yes" or "no".');
end

runCmd = InputFileElement();
runCmd.createInputFileText = @(~) thermalVelocitiesCfg(temperature, zeroTotalMom, seed);

end

