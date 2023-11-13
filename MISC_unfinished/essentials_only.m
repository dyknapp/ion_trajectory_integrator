% Variables to set beforehand

simion_path = "SIM001_data/001_linear_paul_trap";

electrode_names = ["essentials_only_rf.patxt", ...
                   "essentials_only_dc1.patxt", ...
                   "essentials_only_dc2.patxt", ...
                   "essentials_only_dc3.patxt", ...
                   "essentials_only_dc4.patxt", ...
                   "essentials_only_dc5.patxt", ...
                  ];

start_line = 19;

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

d = 1.0;                % mm

start_time =  0.0;      % us
end_time   =  250.0;     % us

m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  0.001;      % mm

%% RF and endcap electrode voltages.

time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);

step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

RF_frequency = 10.0e+6;
RF_amplitude = 100.0;
endcap_voltage = 10.0;

% RF electrodes: 1
voltages(:, 1) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);

% DC electrode pairs: 2-6
voltages(:, 2) = endcap_voltage * ones([1 time_steps]);
voltages(:, 6) = endcap_voltage * ones([1 time_steps]);

%% Initializing

xx1 = d * double(dimensions(1)) / 2.0; % mm
yy1 = d * double(dimensions(2)) / 2.0; % mm
zz1 = d * double(dimensions(3)) / 2.0; % mm

T = 1;
maxwell = @(v) maxwell_pdf(v, m, T);
v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
theta = 2 * pi * rand();
phi = 2 * pi * rand();
vxx1 = v * sin(theta) * cos(phi);      % mm / us
vyy1 = v * sin(theta) * sin(phi);      % mm / us
vzz1 = v * cos(theta);                 % mm / us

% Hydrogen speed at 2.5eV, small variation
cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);
vzz1 = normrnd(1.0e-3 * sqrt(2 * 2.5 * cmr), 0.1);

%% Integration: FORTRAN version
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
    = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);
its = int32(its);
original_length = length(ts);
x_traj = x_traj(1:its);
y_traj = y_traj(1:its);
z_traj = z_traj(1:its);
ts     =     ts(1:its);
exs    =    exs(1:its);
eys    =    eys(1:its);
ezs    =    ezs(1:its);

elapsed_time = toc;

fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(its / elapsed_time));
%% Plot result

stacked_motion_plots(true, ts, x_traj, y_traj, z_traj, ...
    exs, eys, ezs, RF_frequency, RF_amplitude, start_time, elapsed_time)



