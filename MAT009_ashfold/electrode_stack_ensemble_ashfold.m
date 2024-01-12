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

constants = physical_constants();

%% Variables to set beforehand
%  This code block contains all scalar variables to be set before the
%   simulation is started
% IMPORTANT: electrode voltages are set in the NEXT code block

% Directory information
% Path to simion folder containing potential files *.patxt
simion_path = "D:\ownCloud\01_experiment\07_SIM001_data\009_ashfold_vmi_stack";

% Names of individual files containing electrodes' potentials
electrode_names = ["ion_optics_assy-ashfold-A0.0.pa1.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa2.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa3.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa4.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa5.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa6.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa7.patxt", ...
                   "ion_optics_assy-ashfold-A0.0.pa8.patxt", ...
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
end_time   =  1.0;     % us

% Particle specifications
m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

% Integration step size limit
% The integration scheme used in this software has variable timesteps.  The
%   timesteps will be adjusted so that whatever value is set below, the
%   particle should never travel further than that in one timestep.
maxdist =  0.01;      % mm
n_particles = 1e3;  % 2^10;
cloud_radius = 0.5;
T = 0.1;

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
voltages = set_voltage_at_time(1, 25.0, turn_on_time, step_times, voltages);
% Always have the detector at -2.5kV
voltages = set_voltage_at_time(8, -2500.0, 0.0, step_times, voltages);

% If you are curious to see what the potential from electrode1 would be:
% potential_maps_reshaped = reshape(potential_maps, [length(electrode_names) dimensions]);
% imagesc(squeeze(potential_maps_reshaped(2,round(dimensions(1)/2),:,:)))
% axis image
% colorbar()
% colormap('turbo')
% clear potential_maps_reshaped
%% Initializing
% We need to choose the particle's initial phase space position.

% Initialize the particle in the center of the trap.
% We multiply by d (the grid spacing) because the integrator expects
%   physical units, not indices.
xx1 = d * double(dimensions(1) + 1) / 2.0; % mm
yy1 = d * double(dimensions(2) + 1) / 2.0; % mm
zz1 = d * double(dimensions(3) + 1) * 0.945; % mm

% If you are curious to see what the potential from electrode1 would be, 
% compared to the starting point:
potential_maps_reshaped = reshape(potential_maps, [length(electrode_names) dimensions]);

% tiledlayout(4, 2)
% for i = 1:length(electrode_names)
%     nexttile
%     imagesc(squeeze(potential_maps_reshaped(i,round(dimensions(1)/2),:,:)))
%     title(sprintf('Electrode %d', i))
%     axis image
%     colorbar()
%     colormap('turbo')
%     yline(yy1/d, 'r'); xline(zz1/d, 'r');
%     title(sprintf('Electrode %d', i))
% end
% drawnow('update')


vxxs = zeros([1 n_particles]);
vyys = zeros([1 n_particles]);
vzzs = zeros([1 n_particles]);
xxs  = xx1 + normrnd(0, cloud_radius, [1 n_particles]);
yys  = yy1 + normrnd(0, cloud_radius, [1 n_particles]);
zzs  = zz1 + normrnd(0, cloud_radius, [1 n_particles]);

fragment_energy_eV = 0.6;
dissoc_speed = 1.0e-3 * sqrt(2 * fragment_energy_eV * constants("elementary charge") ...
                                / constants("proton mass")); %mm / us
% If we put the ion at the center of the trap with zero initial speed, it
%   will just sit there quietly.  A nice way to give it some speed is by
%   specifying a temperature and choosing a speed based on that.
maxwell = @(v) maxwell_pdf(v, m, T);
% Maxwell-Boltzmann distribution speeds
vs = general_distribution(n_particles, 0.01, 10000, maxwell) * 1.0e-3;
% Uniform distribution direcitions
thetas = 2 * pi * rand([n_particles 1]);
phis   = 2 * pi * rand([n_particles 1]);
% Convert speed & direction -> velocity
vxxs = vs .* sin(thetas) .* cos(phis);      % mm / us
vyys = vs .* sin(thetas) .* sin(phis);      % mm / us
vzzs = vs .* cos(thetas);                   % mm / us

cos2 = @(theta) (cos(theta).^2.) / (pi);
thetas2 = general_distribution(n_particles, 0.001, 2*pi, cos2);
vxxs = vxxs + dissoc_speed * cos(thetas2);
vyys = vyys + dissoc_speed * sin(thetas2);

%% Integration
fprintf("Simulation started.\n")
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[xss, yss, zss, tss, itss] = fly_ensemble(int32(1000), int32(n_particles), xxs, yys, zzs, vxxs, vyys, vzzs, ...
                                      potential_maps, voltages, step_times, ...
                                      int32(time_steps), dimensions, int32(is_electrode), ...
                                      int32(length(electrode_names)), m, q, d, ...
                                      maxdist, end_time);
elapsed_time = toc;
fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(sum(itss) / elapsed_time));

%%

figure;
image_res = 1024;
point_size = 3;
image = zeros(image_res);
coords = zeros([image_res image_res 2]);
for idx1 = 1:image_res
    for idx2 = 1:image_res
        coords(idx1, idx2, :) = d * [double(idx1 * dimensions(1)) / double(image_res), double(idx2 * dimensions(2)) / double(image_res)];
    end
end
parfor idx = 1:n_particles
    image = image + exp((-double(coords(:, :, 1) - xss(idx, end)).^2 ...
                         -double(coords(:, :, 2) - yss(idx, end)).^2) / 0.5);
end
imshow(image / max(image(:)));

% %%
% figure
% [N,C] = hist3([xss(:, end), zss(:, end)],'CDataMode','auto');
% wx = C{1}(:);
% wy = C{2}(:);
% H = pcolor(wx, wy, N');
% set(H,'edgecolor','none');
% view(2)
% shading interp
% axis image
% box on
% colormap jet
% 
% % contourf(C{1}, C{2}, N)

%%


%%

xlims = [-50 50];
ylims = [-10 360];

aspect_ratio_figure = 1.5;
width_screen_pixels = 780;
height_screen_pixels = width_screen_pixels*aspect_ratio_figure;
resolution_dpi = 300;

aspect_ratio_axes = 3.7;
font_size = 12;
line_width_axes = 1;
line_width_plot = 1;


%%% Figure defintion
fh.ion_optics_ashfold = figure(...
    Name                    = 'plot_ion_optics_ashfold',...
    WindowStyle             = 'normal',...
    WindowState             = 'normal',...
    color                   = 'w',...
    Units                   = 'pixels',...
    Position                = [0 0 width_screen_pixels height_screen_pixels],...
    Resize                  = 'off'...
    );

%%%
tl.ion_optics_ashfold = tiledlayout(...
    fh.ion_optics_ashfold,1,1,...
    Padding                 = 'compact',...
    TileSpacing             = 'tight' ...
    );

%%%
ax.ion_optics_ashfold = axes(tl.ion_optics_ashfold);
ax.ion_optics_ashfold.PlotBoxAspectRatio = [1 aspect_ratio_axes 1];
ax.ion_optics_ashfold.Box = 'on';
ax.ion_optics_ashfold.Color = 'none';
ax.ion_optics_ashfold.TickLabelInterpreter = 'latex';
ax.ion_optics_ashfold.FontSize = font_size;
ax.ion_optics_ashfold.XLim = xlims;
ax.ion_optics_ashfold.YLim = ylims;
ax.ion_optics_ashfold.XScale = 'linear';
ax.ion_optics_ashfold.YScale = 'linear';
ax.ion_optics_ashfold.XMinorTick = 'on';
ax.ion_optics_ashfold.YMinorTick = 'on';
ax.ion_optics_ashfold.XGrid = 'off';
ax.ion_optics_ashfold.YGrid = 'off';
ax.ion_optics_ashfold.LineWidth = line_width_axes;
ax.ion_optics_ashfold.YColor = [0 0 0];

xlabel( ax.ion_optics_ashfold,    '$x~$[mm]',                'Interpreter','latex' );
ylabel( ax.ion_optics_ashfold,    '$z~$[mm]',    'Interpreter','latex' );      % '$\alpha~(e^2a_0^2E_H^{-1})$'
hold(   ax.ion_optics_ashfold,    'on' );



%%%

horizontal_offset = -zz1;
vertical_offset = -45.5;

electrodes = ~flipud( squeeze( is_electrode(:, round(dimensions(2)/2.), :) ) )';

imagesc( ...
    ax.ion_optics_ashfold, ...
    d*double(1:dimensions(1))+vertical_offset, ...
    d*double(1:dimensions(3)), ...
    cat( 3, electrodes, electrodes, electrodes ) ...
    )


plot( ...
    ax.ion_optics_ashfold, ...
    yss'+vertical_offset, ...
    zss', ...
    Color = [0 0 0], ...
    LineStyle='-', ...
    LineWidth=0.25 ...
    );


%%

plot3( ...
    xss'+vertical_offset, ...
    yss'+vertical_offset, ...
    zss', ...
    Color = [0 0 0], ...
    LineStyle='-', ...
    LineWidth=0.25 ...
    );


%%


