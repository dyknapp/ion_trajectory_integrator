function [ runCmd ] = evolve( nSteps )
%EVOLVE evolves the simulation for the given number of timesteps.
% 
% Syntax: evolve(nSteps);
% 
% See Also: http://lammps.sandia.gov/doc/run.html

runCmd = InputFileElement();
runCmd.createInputFileText = @(~) runCommandCfg(nSteps);

end