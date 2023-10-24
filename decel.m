%%

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    load("decel_potential.mat");
end
electrode_names = ["einzel3" "einzel2" "einzel1" "endcap1" "endcap2" "rfelecs" "middle1" "middle2" "middle3"];

RF_frequency = 13.2e+6;
RF_amplitude = 270.0;
offset_volts = 0.0;
endcap_volts = 10.0;

m = 2 * 1.00727647 + 0.00054858;
q = 1.0;

end_time = 75.0;
maxdist  = 2.5e-4;

d = 0.125;

xx1 =      d * double(dimensions(1) + 1) / 2.0;
yy1 =      d * double(dimensions(2) + 1) / 2.0;
zz1 =      150.0;

vxx1 = normrnd(0.0, 0.05);
vyy1 = normrnd(0.0, 0.05);
vzz1 =  -10.0;

%%

% RF and endcap electrode voltages.
time_steps_per_us = 100;
time_steps = round((end_time - 0.0) * time_steps_per_us);
step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

rf_electrode = 6;
endcaps_in = 4;
endcaps_out = 5;

ton = 9.0;

voltages(:, rf_electrode) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);
voltages(:, rf_electrode) = voltages(:, rf_electrode) + offset_volts;
voltages(step_times > ton, endcaps_in) = endcap_volts * ones([sum(step_times > ton) 1]);
voltages(:, endcaps_out) = endcap_volts * ones([time_steps 1]);

share = 0.21;
voltages(step_times > ton, 7) = share * endcap_volts * ones([sum(step_times > ton) 1]);
voltages(:, 8) = share * endcap_volts * ones([time_steps 1]);
voltages(:, 9) = share * endcap_volts * ones([time_steps 1]);

voltages(:, 2) = -50.0 * ones([time_steps 1]);

% Simulation
fprintf("RF frequency = %.1f MHz\n" + ...
        "RF amplitude = %.1f V\n" + ...
        "DC offset    = %.1f V\n", ...
        RF_frequency / 1.0e+6, ...
        RF_amplitude, ...
        offset_volts);

potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);

tic
[x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
    = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);
toc
    
its = int32(its);
original_length = length(ts);
x_traj = x_traj(1:its);
y_traj = y_traj(1:its);
z_traj = z_traj(1:its);
ts     =     ts(1:its);
exs    =    exs(1:its);
eys    =    eys(1:its);
ezs    =    ezs(1:its);

%%

particles = 1;

xs = normrnd(0.0, 0.1, [particles, 1]) + (d * double(dimensions(1) + 1) / 2.0);
ys = normrnd(0.0, 0.1, [particles, 1]) + (d * double(dimensions(2) + 1) / 2.0);
zs = normrnd(0.0, 0.1, [particles, 1]) + 150.0;
vxs = normrnd(0.0, 0.1, [particles, 1]);
vys = normrnd(0.0, 0.1, [particles, 1]);
vzs = normrnd(0.0, 0.1, [particles, 1]) - 10;

[x_trajs, y_trajs, z_trajs, tss, exss, eyss, ezss, itss] ...
    = ensemble_trajectory_integration_module(2^20, particles, xs, ys, zs, vxs, vys, vzs, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);

%%

plot(ts, z_traj)
xline(ton, '--k')
yline(25, '-k')
yline(25 + 48, '-k')
legend('Z trajectory', 'Entrance endcap turn-on', 'Endcap positions')
title('Injected H2+ ion, Z trajectory')
ylabel('z (mm)')
xlabel('t (us)')

%%

% plot(d * (1:dimensions(3)), 10.0 * squeeze(potential_maps(5, dimensions(1)/2, dimensions(2)/2, :)))
% plot(d * (1:dimensions(3)), 10.0 * (squeeze(potential_maps(5, dimensions(1)/2, dimensions(2)/2, :))...
%     + 0.0 * squeeze(potential_maps(7, dimensions(1)/2, dimensions(2)/2, :)) ...
%     + 0.25 * squeeze(potential_maps(8, dimensions(1)/2, dimensions(2)/2, :)) ...
%     + 0.25 * squeeze(potential_maps(9, dimensions(1)/2, dimensions(2)/2, :))))
plot(d * (1:dimensions(3)), 10.0 * (squeeze(potential_maps(5, dimensions(1)/2, dimensions(2)/2, :)) + squeeze(potential_maps(4, dimensions(1)/2, dimensions(2)/2, :))...
    + 0.25 * squeeze(potential_maps(7, dimensions(1)/2, dimensions(2)/2, :)) ...
    + 0.25 * squeeze(potential_maps(8, dimensions(1)/2, dimensions(2)/2, :)) ...
    + 0.25 * squeeze(potential_maps(9, dimensions(1)/2, dimensions(2)/2, :))))
xlabel('Position (mm)')
ylabel('Potential (V)')
title('Axial deceleration potential: after capture')