% Basic example illustrating the operation of a hyperbolic paul trap
% r0^2 = 2mm, z0^2 = 1mm
% Trap is a solid of rotation with respect to the x axis
% See ../001_linear_paul_trap/linear_paul_trap.m for substantive comments
%
%  dknapp,  1.9.2023: Wrote script

%% Variables to set beforehand
simion_path = "SIM001_data/002_hyperbolic_paul_trap";
electrode_names = ["hyperbolic_trap.pa1.patxt", ...
                   "hyperbolic_trap.pa2.patxt"
                  ];
start_line = 22;
loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

d = 0.1;                % mm
start_time =  0.0;      % us
end_time   =  50.0;     % us
m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  1.0e-5;    % mm

%% RF and endcap electrode voltages.
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);


% in Hz
RF_frequency = 10.0 * 10.0^6;
% in volts
RF_amplitude = -8.2;
endcap_voltage = 0.0;

% RF electrode: 2
rf_electrode = 2;
voltages(:, rf_electrode) = RF_amplitude * ...
    sin(2 * pi * step_times * RF_frequency / 10.0^6);
voltages(:, rf_electrode) = voltages(:, rf_electrode) + endcap_voltage;

%% Initializing
xx1 = d * double(dimensions(1) - 1) / 2.0; % mm
yy1 = d * double(dimensions(2) - 1) / 2.0; % mm
zz1 = d * double(dimensions(3) - 1) / 2.0; % mm

T = 5;
maxwell = @(v) maxwell_pdf(v, m, T);
v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
theta = 2 * pi * rand();
phi = 2 * pi * rand();
vxx1 = v * sin(theta) * cos(phi);      % mm / us
vyy1 = v * sin(theta) * sin(phi);      % mm / us
vzz1 = v * cos(theta);                 % mm / us

%% Integration: FORTRAN version
tic
[x_traj, y_traj, z_traj, ts, its, recorded] ...
    = integrate_trajectory_lite(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                      potential_maps(:), voltages(:), step_times, ...
                      time_steps, dimensions, int32(is_electrode(:)), ...
                      length(electrode_names), m, q, d, maxdist, end_time, 0.01);

its = int32(its);
recorded = int32(recorded);
original_length = length(ts);
x_traj = x_traj(1:recorded);
y_traj = y_traj(1:recorded);
z_traj = z_traj(1:recorded);
ts     =     ts(1:recorded);

elapsed_time = toc;

fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(its / elapsed_time));



