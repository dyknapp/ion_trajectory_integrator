%% Make a zero potential region for testing

particles = 32;

d = 1.0;
electrode_names = ["zero"];
dimensions = [101 101 101];
potential_maps = reshape(zeros(dimensions), [1 dimensions]);
is_electrode = zeros(dimensions, "int32");

start_time =  0.0;
end_time   =  5.0;

ms = 1.0 * ones([1 particles]);
qs = 1.0 * ones([1 particles]);

maxdist = 1.0e-8;

% Voltages
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

%% Initial conditions
T = 300.;

initial_vector = zeros([1 particles]);
xs  = initial_vector;
ys  = initial_vector;
zs  = initial_vector;
vxs = initial_vector;
vys = initial_vector;
vzs = initial_vector;

% Positions
xs = normrnd(d * double(dimensions(1) - 1) / 2., 10.0, [1 particles]);
xs(xs > d * double(dimensions(1) - 3)) = d * double(dimensions(1) - 3);
xs(xs < d * 3.0)                       = d * 3.;

ys = normrnd(d * double(dimensions(2) - 1) / 2., 10.0, [1 particles]);
ys(ys > d * double(dimensions(2) - 3)) = d * double(dimensions(2) - 3);
ys(ys < d * 3.0)                       = d * 3.;

zs = normrnd(d * double(dimensions(3) - 1) / 2., 10.0, [1 particles]);
zs(zs > d * double(dimensions(3) - 3)) = d * double(dimensions(3) - 3);
zs(zs < d * 3.0)                       = d * 3.;

% Maxwell-Boltzmann
for idx = 1:particles
    maxwell = @(v) maxwell_pdf(v, ms(idx), T);
    v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
    theta = 2 * pi * rand();
    phi = 2 * pi * rand();
    vxs(idx) = v * sin(theta) * cos(phi);      % mm / us
    vys(idx) = v * sin(theta) * sin(phi);      % mm / us
    vzs(idx) = v * cos(theta);                 % mm / us
end

% scatter3(xs, ys, zs, '.');
% axis image;
% xlim([0 d*dimensions(1)])
% ylim([0 d*dimensions(2)])
% zlim([0 d*dimensions(3)])

%% Integration
% Simulation params
record_interval = 1.0e-3;
interps         = int32(1.0e+3);

tic
[x_traj, y_traj, z_traj, ts, its] ...
    = fly_cloud(record_interval,  int32(interps), int32(particles),...
                xs, ys, zs, vxs, vys, vzs, ...
                potential_maps, voltages, step_times, ...
                int32(time_steps), int32(dimensions), int32(is_electrode), ...
                int32(length(electrode_names)), ms, qs, d, maxdist, end_time);
elapsed = toc;
fprintf("Simulation finished ( %5.1f s).  Iterations: %9.1f ( %.3g its/s)\n", elapsed, its, its / elapsed)



%%

plot(ts, x_traj)






