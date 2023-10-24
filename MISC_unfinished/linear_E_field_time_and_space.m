
%  dknapp,  18.9.2023: Wrote script

%% Variables to set beforehand

d = 1.0;                % mm
start_time =  0.0;      % us
end_time   =  100.0;    % us
m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  1.0e-6;    % mm

%% Create the potential.
%  d = 1.0, so 1mm/gu
dimensions = int32([11 11 104]);
potential_maps = zeros([2 dimensions]);
is_electrode = zeros(dimensions);

field_strength = 100.0; % V/m
for z = 1:dimensions(3)
    for x = 1:dimensions(1)
        for y = 1:dimensions(2)
            potential_maps(1, x, y, z) = (1.0) *  0.5 * double(z) * field_strength;
            potential_maps(2, x, y, z) = (1.0) * -0.5 * double(z) * field_strength;
        end
    end
end

%% Electrode voltages.
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, 2]);

idx = 0
for t = step_times
    idx = idx + 1;
    voltages(idx,1) = -1.0 * t;
    voltages(idx,2) =  1.0 * t;
end

%% Initializing
xx1 = d * double(dimensions(1) - 1) / 2.0; % mm
yy1 = d * double(dimensions(2) - 1) / 2.0; % mm
zz1 = 2.0; % mm

vxx1 = 0.01;      % mm / us
vyy1 = 0.02;      % mm / us
vzz1 = 0.0;       % mm / us

%% Integration: FORTRAN version
tic
[x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
    = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      int32(2), m, q, d, maxdist, end_time);

its = int32(its);
original_length = length(ts);
x_traj = x_traj(1:its);
y_traj = y_traj(1:its);
z_traj = z_traj(1:its);
ts     =     ts(1:its);
exs    =    exs(1:its);
eys    =    eys(1:its);
ezs    =    ezs(1:its);

if its == 1048576
    fprintf("!! Max allowed number of timesteps reached.                  !!\n" + ...
            "!! Theory calculation is based on partial distance traveled. !!\n")
end

elapsed_time = toc;
fprintf("Simulation took:                       %.3g s (%d it/s)\n", elapsed_time, round(its / elapsed_time));
fprintf("Mean timestep:                         %.8f us\n", (ts(end) / double(its)));
fprintf("Time of last data point recorded:      %.8f   us\n", ts(end));
fprintf("Theory prediction:                     %.8f   us\n", ...
    1.0e+6 * ((6 * (1.0e-3) * (z_traj(end) - zz1 * d) * 1.660539067e-27 * m) / (q * 1.602176634e-19 * 1000 * 1.0e+6 * field_strength))^(1.0/3.0))


