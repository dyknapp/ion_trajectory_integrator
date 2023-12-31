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
RF_frequency = 100.0 * 10.0^6;
% in volts
RF_amplitude = 10000.0;
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
% elapsed_time = toc;
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

stacked_motion_plots(true, ts, x_traj, y_traj, z_traj, ...
    exs, eys, ezs, RF_frequency, RF_amplitude, start_time, elapsed_time);

%% FFT analysis
if((abs(ts(end) - end_time) < (end_time - start_time) * 1e-3) ...
        || (its == original_length))
    motion_fft(ts, x_traj, RF_frequency, RF_amplitude, 'x');
end



