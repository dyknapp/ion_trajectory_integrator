function SetSimulationDomain( sim, length, width, height )
%SETSIMULATIONDOMAIN Sets the initial size of the simulation 'box'. The box
%may expand larger than this size during the simulation but all atoms
%should be initially placed in the box. Units are SI (m).
%
% SYNTAX: SetSimulationDomain(sim, length, width, height )
%
% Example:
%  SetSimulationDomain(sim, 1e-3, 1e-3, 1e-3)
sim.SimulationBox = struct('length', length, 'width', width, 'height', height);
end

