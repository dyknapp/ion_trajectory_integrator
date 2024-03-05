function [ output ] = defineSimulationBox(  l, w, h, ionSpecies  )
%DEFINESIMULATIONBOX Defines a simulation box.

output = struct('cfgFileHandle', @()simulationBoxCfg( l, w, h, ionSpecies ));

end