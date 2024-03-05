function [ fix ] = harmonicOscillator( f_x, f_y, f_z, group )
%HARMONICOSCILLATOR Adds a simple harmonic oscillator to the given group of
%atoms (all atoms if no group specified).
% f_x, f_y, f_z are spring constants along the three cartesian axes.
%
% SYNTAX: harmonicOscillator( f_x, f_y, f_z )
%         harmonicOscillator( f_x, f_y, f_z, group )

if nargin > 3
    if ~isa(group, 'LAMMPSGroup')
        error('group argument must be of type LAMMPSGroup');
    end
end

if f_x < 0 || f_y < 0 || f_z < 0
    warning('One or more spring constants are less than zero (repulsive)');
end

fix = LAMMPSFix();
%determine limit on timestep.
fix.time = 1/(20*max([f_x, f_y, f_z]));
fix.createInputFileText = @harmonicOscillatorCfg;
    
if nargin > 3
    fix.InputFileArgs = { f_x, f_y, f_z, group };
else
    fix.InputFileArgs = { f_x, f_y, f_z };
end


end

