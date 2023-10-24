loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    load("test_potential.mat");
end
electrode_names = ["einzel3" "einzel2" "einzel1" "endcap1" "endcap2" "rfelecs"];

RF_frequency = 13.2e+6;
RF_amplitude = 270.0;
offset_volts = 0.0;
endcap_volts = 10.0;

m = 1.00727647 * 2 + 0.00054858;
q = 1.0;

end_time = 15.0;
maxdist  = 2.5e-4;

d = 0.125;

xx1 =      d * double(dimensions(1) + 1) / 2.0;
yy1 =      d * double(dimensions(2) + 1) / 2.0;
zz1 =      d * 400.0;

vxx1 = normrnd(0.0, 0.05);
vyy1 = normrnd(0.0, 0.05);
vzz1 =  0.01;

% RF and endcap electrode voltages.
time_steps_per_us = 100;
time_steps = round((end_time - 0.0) * time_steps_per_us);
step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

rf_electrode = 6;
endcaps_in = 4;
endcaps_out = 5;


voltages(:, rf_electrode) = voltages(:, rf_electrode) + offset_volts;
voltages(:, endcaps_in) = endcap_volts * ones([time_steps 1]);
voltages(:, endcaps_out) = endcap_volts * ones([time_steps 1]);
for x = 1:dimensions(1)
    potential_maps(endcaps_in, x, :, :) = potential_maps(endcaps_in, x, :, :) + 0.001*double(x);
end

%%
res = 2048;
tss = zeros([res 1]);
depths = tss;
amps = tss;
tfs = tss;
idx = 0;
rfas = linspace(0, 600, res);
for RF_amplitude = rfas
    idx = idx + 1;
    voltages(:, rf_electrode) = RF_amplitude * ...
        cos(2 * pi * step_times * RF_frequency / 10.0^6);
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
    disp(ts(its));
    tss(idx) = ts(its);
    original_length = length(ts);
    x_traj = x_traj(1:its);
    y_traj = y_traj(1:its);
    z_traj = z_traj(1:its);
    ts     =     ts(1:its);
    exs    =    exs(1:its);
    eys    =    eys(1:its);
    ezs    =    ezs(1:its);
    
    motion_data = x_traj;
    clean_traj = motion_data - mean(motion_data);
    % Resample so that points are equally spaced
    new_samples = linspace(ts(1), ts(end), ...
                           round(length(clean_traj) / 2) * 2);
    dt = new_samples(2) - new_samples(1);
    clean_traj = interp1(ts, clean_traj, new_samples);
    
    L = double(length(clean_traj));
    P2 = abs(fft(clean_traj) / L);
    P1 = P2(1:L / 2 + 1);
    P1(2:end - 1) = 2 * P1(2:end - 1);
    
    f = (1.0 / dt) * (0:(L/2))/L;
    % figure
    % plot(f,P1, '-')
    % xlim([0, RF_frequency * 2.0e-6])
    % xlabel('Frequency (MHz)')
    % ylabel('Amplitude (mm)')
    % title(sprintf('Frequency components of ion motion %s component', motion_name))
    % xline(RF_frequency * 1.0e-6, '--')
    max_index = find(P1 == max(P1), 1);
    amps(idx) = max(P1);
    max_freq = f(max_index) * 1.0e+6;
    % xline(RF_frequency * 1.0e-6 - max_freq * 1.0e-6, '--')
    % xline(RF_frequency * 1.0e-6 + max_freq * 1.0e-6, '--')
    % legend('FFT', 'RF Frequency', 'RF - Secular', 'RF + Secular')
    
    % Find peak frequency component -> calculate "harmonic" potential shape
    fprintf('Secular x frequency:                 %.3g MHz\n', max_freq / 1.0e+6);
    tfs(idx) = max_freq / 1.0e+6;
    
    % r_0 calculation -> see notes (???)
    r0 = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27^2 * RF_frequency^2 * max_freq^2))^0.25;
    fprintf("r_0 predicted from FFT result is:    " + ...
            "%.3g mm\n", 1.0e+3 * r0);
    
    % Trap depth
    trap_depth = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27 * RF_frequency^2 * r0^2)) / 1.602e-19;
    depths(idx) = trap_depth;
    fprintf("Trap depth is:                       " + ...
            "%.3g eV (%.3g K)\n", trap_depth, trap_depth / 8.617e-5);
end

%%
hold on
% plot(rfas, depths);
plot(rfas, amps);

hold off


%%

otss = tss;
odepths = depths;
oamps = amps;
otfs = tfs;