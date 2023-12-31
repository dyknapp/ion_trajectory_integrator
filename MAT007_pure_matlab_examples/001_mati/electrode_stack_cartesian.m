% This is a basic example illustrating a MATI using an electrode stack
% This is done in Cartesian coordinates; for cylindrical, see the
% accompanying example.
%
% A lot of the comments below will be copied verbatim from the example:
%                   MAT002_examples/001_linear_paul_trap/linear_paul_trap.m
% This is the main example for this repository, and it is more likely that
% I remember to update it with new changes, especially when it comes to
% further FORTRAN development.  The MATLAB side will probably stay static.
% The main purpose of reproducing the comments here is for convenience.
%
% 24.10.2023,  dknapp: wrote the example

%% Variables to set beforehand
%  This code block contains all scalar variables to be set before the
%   simulation is started
% IMPORTANT: electrode voltages are set in the NEXT code block

% Directory information
% Path to simion folder containing potential files *.patxt
simion_path = "SIM001_data/001_linear_paul_trap";

% Names of individual files containing electrodes' potentials
electrode_names = ["backingplate.patxt", ...
                   "electrode1.patxt", ...
                   "electrode2.patxt", ...
                   "electrode3.patxt", ...
                   "electrode4.patxt", ...
                   "electrode5.patxt", ...
                   "electrode6.patxt", ...
                   "detector.patxt", ...
                  ];

% Where does the data in the patxt file start?
start_line = 19;

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
    % potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

% What is the grid spacing in the *.patxt files that were loaded?
d = 1.0;                % mm

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
maxdist =  0.001;      % mm

%% Electrode voltages.
% I wrote a simplified system for setting the electrode voltages.  If you
% don't care about the mechanics of how it's done, just skip to the end of
% this section

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

% To set a static voltage at a certain time, use this function:
% set_voltage_at_time(electrode, voltage, time, step_times, voltages_array)
%
% After turn_on_time (us), set the potential of electrode1.patxt to 25V
turn_on_time = 0.1;
voltages = set_voltage_at_time(2, 25.0, turn_on_time, step_times, voltages);
% Always have the detector at -2.5kV
voltages = set_voltage_at_time(8, -2500.0, 0.0, step_times, voltages);

% If you are curious to see what the potential from electrode1 would be:
potential_maps_reshaped = reshape(potential_maps, [length(electrode_names) dimensions]);
imagesc(squeeze(potential_maps_reshaped(2,round(dimensions(1)/2),:,:)))
axis image
colorbar()
colormap('turbo')
clear potential_maps_reshaped
%% Initializing
% We need to choose the particle's initial phase space position.

% Initialize the particle in the center of the trap.
% We multiply by d (the grid spacing) because the integrator expects
%   physical units, not indices.
xx1 = d * double(dimensions(1) + 1) / 2.0; % mm
yy1 = d * double(dimensions(2) + 1) * 0.9; % mm
zz1 = d * double(dimensions(3) + 1) / 2.0; % mm

% If you are curious to see what the potential from electrode1 would be, 
% compared to the starting point:
potential_maps_reshaped = reshape(potential_maps, [length(electrode_names) dimensions]);

tiledlayout(4, 2)
for i = 1:length(electrode_names)
    nexttile
    imagesc(squeeze(potential_maps_reshaped(i,:,:,round(dimensions(3)/2))))
    title(sprintf('Electrode %d', i))
    axis image
    colorbar()
    colormap('turbo')
    yline(zz1/d, 'r'); xline(yy1/d, 'r');
    title(sprintf('Electrode %d', i))
end
drawnow('update')

% If we put the ion at the center of the trap with zero initial speed, it
%   will just sit there quietly.  A nice way to give it some speed is by
%   specifying a temperature and choosing a speed based on that.
T = 0.1;
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
% fprintf("Simulation started.\n")
% tic
% ion_trajectory =  integrate_trajectory(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
%                                       potential_maps, voltages, step_times, ...
%                                       dimensions, is_electrode, m, q, d, ...
%                                       maxdist, end_time);
% x_traj = ion_trajectory.x;
% y_traj = ion_trajectory.y;
% z_traj = ion_trajectory.z;
% ts     = ion_trajectory.t;
% exs    = ion_trajectory.ex;
% eys    = ion_trajectory.ey;
% ezs    = ion_trajectory.ez;
% its    = length(ts);
% elapsed_time = toc;
% fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(its / elapsed_time));

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

%%

figure
hold on
plot(ts, y_traj);
plot(ts, eys);
legend('Y-trajectory', 'Y E-field', 'Location', 'southwest')
ylabel('Trajectory (mm) & E-field (V/m)')
xlabel('Time (us)')
title(sprintf('Trajectory of ion accelerated at %.1gus', turn_on_time))
xline(turn_on_time);
hold off

%%

figure
imagesc((0:dimensions(2)-1)*d, (0:dimensions(3)-1)*d, squeeze(is_electrode(round(dimensions(1)/2), :, :))')
xlabel("y (mm)"); ylabel("z (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image

hold on
plot(y_traj, z_traj, '-r');
hold off
