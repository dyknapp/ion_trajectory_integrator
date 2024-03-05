function [ runCmd ] = minimize( etol, ftol, maxiter, maxeval, maxDist )
%MINIMIZE Perform an energy minimization of the system by iteratively
%adjusting atom coordinates.
% The system is minimized by evolving the equations of motion with a large
% viscous damping, which rapidly removes energy from the system. The
% maximum distance an atom may move in an iteration step is determined by
% maxDist. The minimisation is terminated when one of the following
% criteria is met:
%
% 1. energy difference (as fraction of absolute energy) between steps is less than etol.
% 2. total force is less than ftol
% 3. the number of iterations exceeds maxiter
% 4. the number of force evaluations exceeds maxeval
%
% SYNTAX: minimize( etol, ftol, maxiter, maxeval, maxDist )
%
% See Also: http://lammps.sandia.gov/doc/minimize.html

runCmd = InputFileElement();
runCmd.createInputFileText = @(~) minimizeCfg(etol, ftol, maxiter, maxeval, maxDist);

end