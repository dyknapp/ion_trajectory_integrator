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

% Place protons d apart + impact factor
b = 0;%1.0e-6;
d = 3.0e-5;
seed = 1;
%                                      x              y           z
cloud = placeAtoms(sim, protons, [-d/2. d/2.]', [-b/2. b/2.]', [0. 0.]');

%Configure outputs.
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 1));
pylion_steps = round(1.0e-6 * end_time / sim.TimeStep);
fprintf("Pylion timesteps: %d\n", pylion_steps);

%% Run simulation
sim.Add(evolve(pylion_steps));
tic
sim.Execute();
toc 

%% Load the data
% Load the results from the output file:
[timestep, ~, xs, ys, zs] = readDump('positions.txt');

%% Plot it

plot(sim.TimeStep * double(timestep), xs(2, :), '-');

%% Do the same calculation with my FORTRAN code
% Protons
particles = 2;
ms = 1.0 * ones([1 particles]);
qs = 1.0 * ones([1 particles]);

% Dummy potential
d_grid = 1.0e-2;
electrode_names = ["one"];
dimensions = [101 101 101];
potential_maps = reshape(zeros(dimensions), [1 dimensions]);
is_electrode = zeros(dimensions, "int32");

% Center offset (we can't center around zero like with LAMMPS)
offset = d_grid * [double(dimensions(1) - 1) / 2., double(dimensions(2) - 1) / 2., double(dimensions(3) - 1) / 2.];

% Voltages
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

% Step sizes
maxdist = 1.0e-8;
record_interval = 0.; % This variable does nothing
interps         = int32(1.0e+6);

quick_FORTRAN_to_mex("fly_cloud")
tic
[x_trajs, y_trajs, z_trajs, ts, its] ...
    = fly_cloud(record_interval,  int32(interps), int32(particles) ,...
                (1.0e+3) * [-d/2. d/2.]' + offset(1), ...
                (1.0e+3) * [-b/2. b/2.]' + offset(2), ...
                (1.0e+3) * [0. 0.]' + offset(3), ...
                [0. 0.]',      [0. 0.]',      [0. 0.]', ...
                potential_maps, voltages, step_times, ...
                int32(time_steps), int32(dimensions), int32(is_electrode), ...
                int32(length(electrode_names)), ms, qs, d_grid, maxdist, end_time);
elapsed = toc;
fprintf("Simulation finished ( %5.1f s).  Iterations: %9.1f ( %.3g its/s)\n", elapsed, its, its / elapsed)
% plot(ts, (x_trajs(2,:) - offset(1)));

%% Superimpose plots
figure
hold on
plot(sim.TimeStep * double(timestep), xs(2, :), '-k', LineWidth=5);
plot(ts * 1.0e-6, (x_trajs(2,:) - offset(1)) * 1.0e-3, '-w');
hold off
xlim(1.0e-6 * [min(ts) max(ts)])

%%
syms y(t)
[V] = odeToVectorField((1 * 1.67262192369e-27) * diff(y,2) == ((1.602176634e-19)^2.) / (4. * pi * (8.8541878128e-12) * (2. * y)^2.));
M = matlabFunction(V,'vars', {'t','Y'});
options = odeset('RelTol', 1e-13, 'Stats', 'on', 'MaxStep', 1e-12);
[t, y] = ode89(M, [0 end_time*1.0e-6], [d/2.; 0.], options);
% hold on
% plot(t, y(:, 1), '-');
% hold off

%% Error comparison

pylion_ys = interp1(sim.TimeStep * double(timestep), xs(2, :), t, "linear");
dknapp_ys = interp1(ts * 1.0e-6, (x_trajs(2,:) - offset(1)) * 1.0e-3, t, "linear");

semilogy(t, 100. * abs(pylion_ys - y(:, 1)) ./ abs(y(:, 1)), '-'); hold on
semilogy(t, 100. * abs(dknapp_ys - y(:, 1)) ./ abs(y(:, 1)), '-'); hold off

title("Percentage disagreement with ODE89 Solution")
xlabel("Time (s)")
ylabel("Log scale: Percentage Disagreement (%)")
legend("pylion error", "dknapp error");
