function [ strings ] = runCommandCfg( numberOfSteps )
%RUNCOMMANDCFG Generates config file text for running a LAMMPS simulation.

strings =           { '#Run simulation', sprintf('run %d', numberOfSteps)};
end

