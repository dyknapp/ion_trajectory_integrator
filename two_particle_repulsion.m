%% Make a zero potential region for testing

particles = 2;

d = 1.0;
electrode_names = ["zero"];
dimensions = [101 101 101];
potential_maps = reshape(zeros(dimensions), [1 dimensions]);
is_electrode = zeros(dimensions, "int32");

for i = 1:dimensions(3)
    potential_maps(1, :, :, i) = 0.0e-3 * double(i) * ones(dimensions(1:2));
end

start_time =  0.0;
end_time   =  5.0;

ms = 1.0 * ones([particles 1]);
qs = 1.0 * ones([particles 1]);

maxdist = 1.0e-7;

% Voltages
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = ones([time_steps, length(electrode_names)]);

%% Initial conditions
impact_speed = 0.1;
distance_mm = 0.1;
impact_factor_mm = 1.0e-3;

initial_vector = zeros([1 particles]);
xs  = initial_vector;
ys  = initial_vector;
zs  = initial_vector;
vxs = initial_vector;
vys = initial_vector;
vzs = initial_vector;

% Positions
xs(1) = (d * double(dimensions(1) - 1) / 2.) + (distance_mm / 2.);
xs(2) = (d * double(dimensions(1) - 1) / 2.) - (distance_mm / 2.);
ys(1) = (d * double(dimensions(2) - 1) / 2.) + (impact_factor_mm / 2.);
ys(2) = (d * double(dimensions(2) - 1) / 2.) - (impact_factor_mm / 2.);
zs = (d * double(dimensions(3) - 1) / 2.) * ones(size(zs));

% Speeds
vxs(1) = -impact_speed / 2.;
vxs(2) =  impact_speed / 2.;

constants = physical_constants();
expected_acceleration = qs(1) * qs(2) * (constants("elementary charge")^2.) ...
    / (constants("proton mass") * 4 * pi * constants("epsilon0") * (distance_mm * 1.0e-3)^2.0);
fprintf("Expected initial accelerationdue to Coulomb repulsion: %.15g m/s^2\n", expected_acceleration)

%% Integration
% Simulation params
record_interval = 1.0e-3;
interps         = int32(1.0e+3);

tic
[x_traj, y_traj, z_traj, ts, its] ...
    = fly_cloud(int32(record_interval),  int32(interps), int32(particles),...
                xs, ys, zs, vxs, vys, vzs, ...
                potential_maps, voltages, step_times, ...
                int32(time_steps), int32(dimensions), int32(is_electrode), ...
                int32(length(electrode_names)), ms, qs, d, maxdist, end_time);
elapsed = toc;
fprintf("Simulation finished ( %5.1f s).  Iterations: %9.1f ( %.3g its/s)\n", elapsed, its, its / elapsed)

%%

plot(ts, x_traj)

%%

plot(ts, y_traj)





