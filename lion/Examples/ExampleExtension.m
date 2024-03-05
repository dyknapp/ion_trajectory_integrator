%% Example extension
% Runs a simulation with damping proportional to v^2.
close all
clearvars
sim = LAMMPSSimulation();
sim.SetSimulationDomain(1e-3,1e-3,1e-3);


charge = 0;
mass = 40;
calciumIons = sim.AddAtomType(charge, mass);

N = 100;
seed = 1;
cloud = createIonCloud(sim, 1e-3, calciumIons, N, seed);

% Give initial velocities to the atoms
initialT = 1e-3;
sim.Add(thermalVelocities(initialT, 'no'));

% Add v^2 damping
sim.Add(velocitySquaredDamping(5e4, calciumIons.Group))

% Configure outputs.
sim.Add(dump('velocity.txt', {'id', 'vx', 'vy', 'vz' }, 10));

% Run simulation
sim.Add(evolve(15000));
sim.TimeStep = 1e-8;
sim.Execute();

%% Load the data
% Load the results from the output file:
[t_v, ~, vx,vy,vz] = readDump('velocity.txt');

%% Plot
% Plot the velocities of each atom over time

set(gcf, 'Color', 'w')
semilogy(t_v, (vx.^2+vy.^2+vz.^2).^0.5)
xlabel('timestep')
ylabel('velocity')
