function [ steps ] = getSteps( time, unit )
%GETSTEPS This function returns the number of steps required in the current
% experimental timestep setup to fulfill a given unit of time.
%
% Example: getSteps(100, 'ms'), getSteps(1e-2, 's')
%
% Note: This function has variables that survive beyond the function scope,
% and as such, may be vulnerable to race conditions. It is preferable to
% use time arguments to functions instead. The intended usage is that
% instead of specifying a number of timesteps in the precalculated cfg
% text, a function handle is given instead. Then, when writing the input
% file (which is written sequentially) these function handles are invoked
% for those fixes/runs that require them, and so use the correct number of
% timesteps for the given simulation timestep size. This allows the
% simulation timestep size to vary according to what the current most
% constraining fix is.
persistent timestep;

%To set the timestep, use getSteps(timestep, 'set'); To clear, use
%getSteps(0, 'set');

if strcmp('unit', 'set')
    if time > 0
        timestep = time;
    else
        clear timestep;
    end
else
    
    if ~exist('timestep', 'var')
        error('"timestep" not set. This function should only be invoked during the creation of the config file for LAMMPS, as it is only during this time that the timestep can be determined. Did you mean to pass the function handle?');
    end
    
    if ~exist('unit', 'var')
        unit = 'steps';
    end
    
    switch unit
        case 'steps'
            steps = time;
        case 's'
            steps = time / timestep;
        case 'ms'
            steps = time * 1e-3 / timestep;
        case 'us'
            steps = time * 1e-6 / timestep;
        case 'ns'
            steps = time * 1e-9 / timestep;
    end
    
end
end

