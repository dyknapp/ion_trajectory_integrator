%% Simple example
% This simple example script creates a single species ion trap. 3 ions are
% confined within this trap, and are cooled by interaction with a langevin
% bath. The positions of atoms are output every 10 timesteps, and a secular
% velocity is also output by averaging atomic velocities over each RF
% cycle.
close all
clearvars

%Create an empty experiment.
sim = LAMMPSSimulation();

% Add a simulation box. This determines the region that will be simulated.
% The simulation box may expand beyond these limits, but these specify the
% minimum dimensions of the simulation box.
sim.SetSimulationDomain(1e-3,1e-3,1e-3);

% Add a new atom type:
charge = 1;
mass = 40;
calciumIons = sim.AddAtomType(charge, mass);

% Create a cloud of ions based on this atom type. Seed is used for random
% number generation in placing the atoms.
N = 100;
seed = 1;
cloud = createIonCloud(sim, 1e-3, calciumIons, N, seed);

%Add the linear Paul trap electric field.
RF = 3.85e6;
z0 = 5.5e-3 / 2;
r0 = 7e-3 / 2;
geometricKappa = 0.244;
U0 = 7;
V0 = 300;
sim.Add(linearPT(V0, U0, z0, r0, geometricKappa, RF));

%Add some damping bath
T = 1e-3;
sim.Add(langevinBath(T, 10e-6, sim.Group(calciumIons)));

%Configure outputs.
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 1));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));

% Run simulation
sim.Add(evolve(15000));
sim.Execute();

%% Load the data
% Load the results from the output file:
[timestep, ~, x,y,z] = readDump('positions.txt');

%% Plot
% Plot the final positions of atoms

color = [112 146 190]/255;
clf;
depthPlot(x(:,end), y(:,end), z(:,end), color, [ 20 70 ]);
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;



