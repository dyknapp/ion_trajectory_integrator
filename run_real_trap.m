loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    load("test_potential.mat");
end
electrode_names = ["einzel3" "einzel2" "einzel1" "endcap1" "endcap2" "rfelecs"];

RF_frequency = 13.2e+6;
RF_amplitude = 39.8632;%270.0;
offset_volts = 0.0;
endcap_volts = 10.0;

m = 1.00727647 + 2.01410178 + 0.00054858;
q = 1.0;

end_time = 1000.0;
maxdist  = 2.5e-4;

d = 0.125;

xx1 =      d * double(dimensions(1) + 1) / 2.0;
yy1 =      d * double(dimensions(2) + 1) / 2.0;
zz1 =      d * 400.0;

vxx1 = normrnd(0.0, 0.05);
vyy1 = normrnd(0.0, 0.05);
vzz1 =  15.0;

% RF and endcap electrode voltages.
time_steps_per_us = 100;
time_steps = round((end_time - 0.0) * time_steps_per_us);
step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

rf_electrode = 6;
endcaps_in = 4;
endcaps_out = 5;

voltages(:, rf_electrode) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);
voltages(:, rf_electrode) = voltages(:, rf_electrode) + offset_volts;
voltages(:, endcaps_in) = endcap_volts * ones([time_steps 1]);
voltages(:, endcaps_out) = endcap_volts * ones([time_steps 1]);

% Simulation
fprintf("RF frequency = %.1f MHz\n" + ...
        "RF amplitude = %.1f V\n" + ...
        "DC offset    = %.1f V\n", ...
        RF_frequency / 1.0e+6, ...
        RF_amplitude, ...
        offset_volts);

repeats = 1;
x_trajs = zeros([repeats*1048576 1]);
y_trajs = zeros([repeats*1048576 1]);
z_trajs = zeros([repeats*1048576 1]);
tss     = zeros([repeats*1048576 1]);
next_write = 1;

potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);

for rep = 1:repeats
    disp(rep)

    tic
    [x_traj, y_traj, z_traj, ts, ~, ~, ~, its] ...
        = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          length(electrode_names), m, q, d, maxdist, end_time);
    toc
    
    its = int32(its);
    x_trajs(next_write:next_write + its - 1) = x_traj(1:its);
    y_trajs(next_write:next_write + its - 1) = y_traj(1:its);
    z_trajs(next_write:next_write + its - 1) = z_traj(1:its);
    if rep > 1
        tss(next_write:next_write + its - 1)     = ts(1:its) + tss(next_write - 1);
    else
        tss(next_write:next_write + its - 1)     = ts(1:its);
    end

    next_write = next_write + its;
    if its < 1048576
        break
    end

    xx1 = x_trajs(next_write - 1);
    yy1 = y_trajs(next_write - 1);
    zz1 = z_trajs(next_write - 1);
    vxx1 = (x_trajs(next_write - 1) - x_trajs(next_write - 2)) / (tss(next_write - 1) - tss(next_write - 2));
    vyy1 = (y_trajs(next_write - 1) - y_trajs(next_write - 2)) / (tss(next_write - 1) - tss(next_write - 2));
    vzz1 = (z_trajs(next_write - 1) - z_trajs(next_write - 2)) / (tss(next_write - 1) - tss(next_write - 2));

    voltages(:, rf_electrode) = RF_amplitude * ...
    cos(2 * pi * (step_times + tss(next_write - 1)) * RF_frequency / 10.0^6);
    voltages(:, rf_electrode) = voltages(:, rf_electrode) + offset_volts;
    voltages(:, endcaps_in) = endcap_volts * ones([time_steps 1]);
    voltages(:, endcaps_out) = endcap_volts * ones([time_steps 1]);
end

%%

figure
plot(tss, x_trajs)
ylabel('x (mm)')
xlabel('t (us)')
xlim([0 5000])
