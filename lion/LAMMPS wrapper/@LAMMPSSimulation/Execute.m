function Execute(sim)
% EXECUTE Creates an input file from the configured simulation, launches
% lammps and runs the simulation. Any changes to the simulation after this
% command (eg adding new atom types) will not affect the results of the
% simulation.
%
% SYNTAX: Execute(sim)

%write the input file.
WriteInputFile(sim);

%invoke LAMMPS
runLAMMPS(sim.ConfigFileName);

sim.HasExecuted = true;

end