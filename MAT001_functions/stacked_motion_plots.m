function t = stacked_motion_plots(new_figure, ts, x_traj, y_traj, z_traj, ...
    exs, eys, ezs, RF_frequency, RF_amplitude, start_time, elapsed_time)

if new_figure
    figure('units','normalized','outerposition',[0 0 1 1])
end
t = tiledlayout(7, 1);
nexttile
plot(ts, x_traj);
xlim([start_time, ts(end)])
% ylim([0, dimensions(1) * d])
ylabel('x (mm)')
title(sprintf("%d total timesteps simulated, %.3g(s) taken (%.0f(it/s))", ...
    length(ts), elapsed_time, ...
    double(length(ts)) / elapsed_time))

nexttile
plot(ts, y_traj);
xlim([start_time, ts(end)])
% ylim([0, dimensions(2) * d])
ylabel('y (mm)')

nexttile
plot(ts, z_traj);
xlim([start_time, ts(end)])
% ylim([0, dimensions(3) * d])
ylabel('z (mm)')

nexttile
plot(ts, RF_amplitude * cos(2 * pi * ts * RF_frequency / 10.0^6));
xlim([start_time, ts(end)])
ylabel('RF Voltage (V)')

% E fields
nexttile
plot(ts, exs);
xlim([start_time, ts(end)])
ylabel('Ex (V/m)')

nexttile
plot(ts, eys);
xlim([start_time, ts(end)])
ylabel('Ey (V/m)')

nexttile
plot(ts, ezs);
xlim([start_time, ts(end)])
ylabel('Ez (V/m)')

xlabel('t (us)')

fprintf("Time of last data point recorded:      %.3g us\n", ts(end));

end