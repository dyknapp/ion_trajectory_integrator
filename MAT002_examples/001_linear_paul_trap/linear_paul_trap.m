% Basic example illustrating the operation of a linear paul trap.
% This particular linear paul trap has 8 electrodes.
% SIMION is used to export a separate .patxt file for each electrode.
% Within the software, the potential within the trap is calculated at each
%   given moment by taking linear sums of the contributions from each
%   electrode.  Following the Broerse convention, the default voltage for
%   each electrode as it is saved in the .patxt files is 10kV.
% You will find the .patxt files in ../../SIM001_data/001_linear_paul_trap
% The geometry of the trap is defined by the quad_v1.gem file, in the same
%   directory
%
%
%  dknapp, 30.8.2023: Wrote script.
%  dknapp,  1.9.2023: Added Maxwell-Boltzmann distribution initialization

%% Variables to set beforehand
%  This code block contains all scalar variables to be set before the
%   simulation is started
% IMPORTANT: electrode voltages are set in the NEXT code block

% Directory information
% Path to simion folder containing potential files *.patxt
simion_path = "SIM001_data/001_linear_paul_trap";

% Names of individual files containing electrodes' potentials
electrode_names = ["quad_doubled.PA1.patxt", ...
                   "quad_doubled.PA2.patxt", ...
                   "quad_doubled.PA3.patxt", ...
                   "quad_doubled.PA4.patxt", ...
                   "quad_doubled.PA5.patxt", ...
                   "quad_doubled.PA6.patxt", ...
                   "quad_doubled.PA7.patxt", ...
                   "quad_doubled.PA8.patxt", ...
                  ];

% How many lines to skip at the head of a *.patxt file?
start_line = 22;

% Load and parse the relevant data that is specified above:
%   potential_maps      :       4D array.  Size: 
%                               [# of electrodes, x grid cells, y ", z "]
%                               Contains the potential from each electrode
%                               being individually set to 1V.  The SIMION
%                               data has been loaded and normalized.
%
%   is_electrode        :       Logical 3D aray.  One cell per grid site
%                               specifies whether it is occupied by an
%                               electrode.  Used to determine if an ion has
%                               collided with an electrode
%
%   dimensions          :       1D array.  Elements 1, 2, 3 respectively
%                               represent the number of grid sits along the
%                               x-, y-, and z-directions.
%
% This is typically the slowest part of the code, so it is skipped if there
%   is already a variable "dimensions" in the current workspace.  This can
%   cause errors if you run this code immediately after having worked on a
%   different set of SIMION potentials, so make sure to clear the workspace
%   in such cases.  If you are nervous and willing to wait each time for
%   the files to load, better to set "loadanyways = true;".
loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

% What is the grid spacing in the *.patxt files that were loaded?
d = 0.1;                % mm

% Time settings in microseconds
% If an ion leaves the simulation domain or hits an electrode, it may not
%   survive all the way until end_time
start_time =  0.0;      % us
end_time   =  5.0;     % us

% Particle specifications
m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

% Integration step size limit
% The integration scheme used in this software has variable timesteps.  The
%   timesteps will be adjusted so that whatever value is set below, the
%   particle should never travel further than that in one timestep.
maxdist =  0.000025;      % mm

%% RF and endcap electrode voltages.
% Here, the voltages for the electrodes are specified.
% In order to accomodate arbitrarily time-varying electrode voltages, you
%   are required to provide the software EQUALLY SPACED time steps, at each
%   step of which, the electrode voltages are specified.  These are passed
%   to the integrator in the form of a [# of timesteps, # of electrodes]
%   size matrix.  During the actual integration, the integrator will
%   perform linear interpolation to find the electrode voltage at an
%   arbitrary time.  Keep in mind that the integrator needs to know the
%   voltages at an arbitrary time because the variable timesteps make it
%   unpredictable at which points in time we actually need to know the
%   voltages.  The number of samples chosen really does little to affect
%   performance, so sample abundantly.  Make sure that your time samples
%   are much shorter than whatever transient voltages you might apply.
time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);
% Just to avoid having to check for edge cases in the interpolation code,
%   the actual voltage specification should go beyond the end of end_time,
%   just a bit.
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

% Below here, write whatever code you need for generating the voltages
% in Hz
RF_frequency = 15.0 * 10.0^6;
% in volts
RF_amplitude = 350.0;
endcap_voltage = 10.0;

% RF electrodes: 1, 2
voltages(:, 1) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);
    
voltages(:, 2) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);

% Center electrodes: 4, 7 (pi phase shift by doing cos -> sin)
% voltages(:, 4) = RF_amplitude * ...
%     sin(2 * pi * step_times * RF_frequency / 10.0^6);
% 
% voltages(:, 7) = RF_amplitude * ...
%     sin(2 * pi * step_times * RF_frequency / 10.0^6);

% Endcaps: 3, 5, 6, 8
for electrode = [3 5 6 8]
    voltages(:, electrode) = ones([time_steps, 1], "double") * endcap_voltage;
end

%% Initializing
% We need to choose the particle's initial phase space position.

% Initialize the particle in the center of the trap.
% We multiply by d (the grid spacing) because the integrator expects
%   physical units, not indices.
xx1 = d * double(dimensions(1)) / 2.0; % mm
yy1 = d * double(dimensions(2)) / 2.0; % mm
zz1 = d * double(dimensions(3)) / 2.0; % mm

% If we put the ion at the center of the trap with zero initial speed, it
%   will just sit there quietly.  A nice way to give it some speed is by
%   specifying a temperature and choosing a speed based on that.
T = 1;
maxwell = @(v) maxwell_pdf(v, m, T);
% Maxwell-Boltzmann distribution speeds
v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
% Uniform distribution directions
theta = 2 * pi * rand();
phi = 2 * pi * rand();
% Convert speed & direction -> velocity
vxx1 = v * sin(theta) * cos(phi);      % mm / us
vyy1 = v * sin(theta) * sin(phi);      % mm / us
vzz1 = v * cos(theta);                 % mm / us

%% Integration: MATLAB version
% ion_trajectory =  integrate_trajectory(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
%                                       potential_maps, voltages, step_times, ...
%                                       dimensions, is_electrode, m, q, d, ...
%                                       maxdist, end_time);
%% Integration: FORTRAN version
% Currently limited to a set number of integration steps, specified as
% MAX_TRAJECTORY_POINTS in the FORTRAN and C files.
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
    = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);

% For simplicity, the FORTRAN integrator assumes that it has carried out
% the maximum permissible number of integration steps.  If the particle has
% died or the simulation has ended before this, we need to trim down the
% output to remove the empty sections at the ends.
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

% Hard to tell what's going on with this one
% figure
% plot3(x_traj, y_traj, z_traj)
% axis equal
% xlim([0, dimensions(1) * d])
% ylim([0, dimensions(2) * d])
% zlim([0, dimensions(3) * d])

% figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(7, 1)
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
plot(ts, cos(2 * pi * ts * RF_frequency / 10.0^6));
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

%% FFT analysis
if((abs(ts(end) - end_time) < (end_time - start_time) * 1e-3) ...
        || (its == original_length))
    % Remove constant offset
    clean_traj = x_traj - mean(x_traj);
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
    figure
    plot(f,P1, '-')
    xlim([0, RF_frequency * 2.0e-6])
    xlabel('Frequency (MHz)')
    ylabel('Amplitude (mm)')
    title('Frequency components of ion motion x component')
    xline(RF_frequency * 1.0e-6, '--')
    max_index = find(P1 == max(P1), 1);
    max_freq = f(max_index) * 1.0e+6;
    xline(RF_frequency * 1.0e-6 - max_freq * 1.0e-6, '--')
    xline(RF_frequency * 1.0e-6 + max_freq * 1.0e-6, '--')
    legend('FFT', 'RF Frequency', 'RF - Secular', 'RF + Secular')
    
    % Find peak frequency component -> calculate "harmonic" potential shape
    fprintf('Secular x frequency:                 %.3g MHz\n', max_freq / 1.0e+6);
    
    % r_0 calculation -> see notes (???)
    r0 = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27^2 * RF_frequency^2 * max_freq^2))^0.25;
    fprintf("r_0 predicted from FFT result is:    " + ...
            "%.3g mm\n", 1.0e+3 * r0);
    
    % Trap depth
    trap_depth = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27 * RF_frequency^2 * r0^2)) / 1.602e-19;
    fprintf("Trap depth is:                       " + ...
            "%.3g eV (%.3g K)\n", trap_depth, trap_depth / 8.617e-5);
end

%% Sanity check FFT of Ex
if((abs(ts(end) - end_time) < (end_time - start_time) * 1e-3) ...
        || (its == original_length))
    new_samples = linspace(ts(1), ts(end), ...
                           round(length(clean_traj) / 2) * 2);
    dt = new_samples(2) - new_samples(1);
    Ex = interp1(ts, exs, new_samples);
    L = double(length(Ex));
    P2 = abs(fft(Ex) / L);
    P1 = P2(1:L / 2 + 1);
    P1(2:end - 1) = 2 * P1(2:end - 1);
    figure
    plot(f,P1, '-')
    xlim([0, RF_frequency * 2.0e-6])
    xlabel('Frequency (MHz)')
    ylabel('Amplitude (V/m)')
    title('Frequency components of ion''s felt E field, x component')
    xline(RF_frequency * 1.0e-6, '--')
    xline(RF_frequency * 1.0e-6 - max_freq * 1.0e-6, '--')
    xline(RF_frequency * 1.0e-6 + max_freq * 1.0e-6, '--')
    legend('FFT', 'RF Frequency', 'RF - Secular', 'RF + Secular')
end



