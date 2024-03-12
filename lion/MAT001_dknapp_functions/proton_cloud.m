close all
clearvars
constants = physical_constants();

% Simulation time
end_time   =  10.0;

% Simulation timescales (us)
record_step =  0.001;
burst_time =   0.000;

% Simulation step scale (mm)
maxdist = 1.0e-8;

%Create an empty experiment.
sim = LAMMPSSimulation();

sim.SetSimulationDomain(1e-3,1e-3,1e-3);
sim.TimeStep = 1.0e-9;

% Add a new atom type:
charge = 1;
mass = 1;
protons = sim.AddAtomType(charge, mass);

% Place protons in a gaussian distribution cloud.
particles = 4;
positions  = normrnd(0.0, 1.0e-6, [particles, 3]);
velocities = zeros([particles, 3]);

ms = ones([1 particles]);
qs = ones([1 particles]);

cloud = placeAtoms(sim, protons, ...
                    squeeze(positions(:, 1)), ...
                    squeeze(positions(:, 2)), ...
                    squeeze(positions(:, 3)));

%Configure outputs.
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z', 'vx', 'vy', 'vz'}, 1));
pylion_steps = round(1.0e-6 * end_time / sim.TimeStep);
fprintf("Pylion timesteps: %d\n", pylion_steps);

%% Run simulation
sim.Add(evolve(pylion_steps));
tic
sim.Execute();
toc 
[timestep, ~, xs, ys, zs, vxs, vys, vzs] = readDump('positions.txt');

%%
% No ion trap potential
omega = 1.0;
depth = 0.0;
R = 1.0;


mh = mexhost;
tic;
[trajectories, times, its, recorded] = ...
    feval(mh, 'nbody', ...
        particles, ...
        positions * 1.0e+3, ...
        velocities * 1.0e-3,  ...
        ms,  ...
        qs, ...
        omega,  ...
        depth,  ...
        R,  ...
        end_time, ...
        maxdist, ...
        record_step, ...
        burst_time ...
    );
elapsed = toc;
fprintf("Simulation finished ( %5.1f s).  Iterations: %9.1f ( %.3g its/s). Mean timestep: %.3g us, Recorded points: %d\n", elapsed, its, its / elapsed,  (times(recorded) / double(its)), recorded)

%% Momenta

momenta = vecnorm(squeeze(sum(trajectories(1:recorded, :, 4:6), 2)), 2, 2);
loglog(times(1:recorded), abs(momenta - momenta(2)) / momenta(2))
set(gca,'TickDir','out');  grid("on");


