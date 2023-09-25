am = 0.0;
qm = 0.25;
RF_frequency = 10.0e+6;
particles = 32;
m = 1.0;
q = 1.0;
end_time = 10.0;
maxdist = 1.0e-5;

% a and q parameters
r0 = 1.0e-3;
e  = 1.602176634e-19;

omega = RF_frequency * (2 * pi);
endcap_voltage = 0.5 * am * (((1.660539067e-27) * m) * omega^2 * r0^2) / e;
RF_amplitude   =       qm * (((1.660539067e-27) * m) * omega^2 * r0^2) / e;

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

d = 0.1;
xs = normrnd(0.0, 0.01, [particles, 1]) + (d * double(dimensions(1) + 1) / 2.0);
ys = normrnd(0.0, 0.01, [particles, 1]) + (d * double(dimensions(2) + 1) / 2.0);
zs = normrnd(0.0, 0.01, [particles, 1]) + (d * double(dimensions(3) + 1) / 2.0);

vxs = normrnd(0.0, 0.001, [particles, 1]);
vys = normrnd(0.0, 0.001, [particles, 1]);
vzs = normrnd(0.0, 0.001, [particles, 1]);

% RF and endcap electrode voltages.
time_steps_per_us = 1000;
time_steps = round((end_time - 0.0) * time_steps_per_us);
step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

% RF electrode: 2
rf_electrode = 2;
voltages(:, rf_electrode) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);
voltages(:, rf_electrode) = voltages(:, rf_electrode) + endcap_voltage;

% Simulation
fprintf("RF frequency = %.1f MHz\n" + ...
        "RF amplitude = %.1f V\n" + ...
        "DC offset    = %.1f V\n" + ...
        "a_r          = %.3g\n" + ...
        "q_r          = %.3g\n", ...
        RF_frequency / 1.0e+6, ...
        RF_amplitude, ...
        endcap_voltage, ...
        am, qm);
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[x_trajs, y_trajs, z_trajs, tss, exss, eyss, ezss, itss] ...
    = ensemble_trajectory_integration_module(particles, xs, ys, zs, vxs, vys, vzs, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);
toc