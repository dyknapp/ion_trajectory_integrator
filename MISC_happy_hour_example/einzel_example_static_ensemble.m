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

constants = physical_constants();

%  This code block contains all scalar variables to be set before the
%   simulation is started
% IMPORTANT: electrode voltages are set in the NEXT code block

% Directory information
% Path to simion folder containing potential files *.patxt
simion_path = "MISC_happy_hour_example";

% Names of individual files containing electrodes' potentials
electrode_names = ["einzel_assembly_shielded.pa1.patxt", ...
                   "einzel_assembly_shielded.pa2.patxt", ...
                   "einzel_assembly_shielded.pa3.patxt", ...
                   "einzel_assembly_shielded.pa4.patxt", ...
                  ];

% Where does the data in the patxt file start?
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
    % potential_maps = potential_maps / 10000.0;
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
maxdist =  0.0001;      % mm

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
time_steps_per_us = 2;
time_steps = round((end_time - start_time) * time_steps_per_us);
% Just to avoid having to check for edge cases in the interpolation code,
%   the actual voltage specification should go beyond the end of end_time,
%   just a bit.
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

voltages(:, 2) = ones([time_steps 1]) * -140.0;

%% Initializing
% We need to choose the particle's initial phase space position.

n_particles = 16;
yy1s = linspace(d * double(dimensions(2) - 1) * 0.4, ...
                d * double(dimensions(2) - 1) * 0.6, ...
                n_particles)';

xx1s  = zeros([n_particles 1]);
zz1s  = zeros([n_particles 1]);
vxx1s = zeros([n_particles 1]);
vyy1s = zeros([n_particles 1]);
vzz1s = zeros([n_particles 1]);

for idx = 1:n_particles
    xx1s(idx) = normrnd(d * double(dimensions(1) + 1) * 0.5, 0.01); % mm
    zz1s(idx) = max(3.01 * d, 6.01 * d + normrnd(0.0, 1.0)); % mm
    
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
    vxx1s(idx) = v * sin(theta) * cos(phi);      % mm / us
    vyy1s(idx) = v * sin(theta) * sin(phi);      % mm / us
    vzz1s(idx) = v * cos(theta);                 % mm / us
    
    % Initial speed
    vz_eV = 1.0;
    vzz1s(idx) = vzz1s(idx) + sqrt(2 * vz_eV * constants("elementary charge") / (m * constants("atomic mass unit")));
end

%% Integration of ensemble in FORTRAN
fprintf("Simulation started.\n")
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[x_trajs, y_trajs, z_trajs, tss, exss, eyss, ezss, itss] ...
    = ensemble_trajectory_integration_module(...
                      2^20, n_particles, xx1s, yy1s, zz1s, vxx1s, vyy1s, vzz1s, ...
                      potential_maps, voltages, step_times, ...
                      time_steps, dimensions, int32(is_electrode), ...
                      length(electrode_names), m, q, d, maxdist, end_time);
elapsed_time = toc;
fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(sum(itss) / elapsed_time));

%%

figure
hold on
for idx = 1:n_particles
    plot(tss(idx, 1:itss(idx)), z_trajs(idx, 1:itss(idx)));
end
hold off
legend('Z-trajectory')
ylabel('Trajectory (mm) & E-field (V/m)')
xlabel('Time (us)')

%%

figure
imagesc(d*(0:dimensions(3)-1), d*(1:dimensions(2)-1), squeeze(is_electrode(:, round(dimensions(2)/2), :)))
xlabel("r (mm)"); ylabel("z (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image
hold on
for idx = 1:n_particles
    plot(z_trajs(idx, 1:itss(idx)), y_trajs(idx, 1:itss(idx)), '-w');
end
hold off

