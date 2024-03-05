close all
clearvars
constants = physical_constants();

% Simulation time
start_time =  0.0;
end_time   =  0.1;

%Create an empty experiment.
sim = LAMMPSSimulation();

sim.SetSimulationDomain(1e-3,1e-3,1e-3);
sim.TimeStep = 1.0e-13;

% Add a new atom type:
charge = 1;
mass = 1;
protons = sim.AddAtomType(charge, mass);

% Place protons in a gaussian distribution cloud.