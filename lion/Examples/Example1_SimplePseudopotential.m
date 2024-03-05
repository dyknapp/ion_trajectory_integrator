close all
clearvars

%Create an empty experiment.
sim = LAMMPSSimulation();

% Add a simulation box. This determines the region that will be simulated.
% The simulation box may expand beyond these limits, but these specify the
% minimum dimensions of the simulation box.
SetSimulationDomain(sim, 1e-3,1e-3,1e-3);

% Define an atom type and add 3 ions of this type to the simulation,
% manually specifying their coordinates:
charge = -1;
mass = 30;
ions = AddAtomType(sim, charge, mass);
placeAtoms(sim, ions, [0 0 0]', [0 0 0]', [-1 0 1]'*1e-4);

%Add the linear Paul trap pseudopotential.
RF = 3.85e6;
sim.Add(linearPseudoPT(300, -0.01, 5.5e-3, 7e-3, 0.244, RF, ions));

%Add a bath to cool the motion of the ions, like for a buffer gas
sim.Add(langevinBath(1e-6, 1e-5));

%Configure outputs. We output position every 10 steps, and also output
%the time-averaged secular velocity every RF cycle.
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 1));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));

% Run simulation
sim.Add(evolve(10000));
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
