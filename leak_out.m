% Variables to set beforehand
simion_path = "";
electrode_names = ["essentials_only_rf.patxt", ...
                   "essentials_only_dc1.patxt", ...
                   "essentials_only_dc2.patxt", ...
                   "essentials_only_dc3.patxt", ...
                   "essentials_only_dc4.patxt", ...
                   "essentials_only_dc5.patxt", ...
                  ];

start_line = 19;

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

d = 1.0;                % mm

start_time =  0.0;      % us
end_time   =  250.0;    % us

m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  0.001;       % mm

n_particles = 1024;
cloud_radius = 0.1;
cloud_length = 3.0;
T = 0.1;
axial_eV = 1.7;

%% RF and endcap electrode voltages.

time_steps_per_us = 1000;
time_steps = round((end_time - start_time) * time_steps_per_us);

step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, length(electrode_names)]);

RF_frequency = 10.0e+6;
RF_amplitude = 250.0;
left_endcap_voltage =  2.5  * 2;
right_endcap_voltage = 1.65 * 2;

% RF electrodes: 1
voltages(:, 1) = RF_amplitude * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);

% DC electrode pairs: 2-6
voltages(:, 2) = left_endcap_voltage * ones([1 time_steps]);
voltages(:, 6) = right_endcap_voltage * ones([1 time_steps]);

%% Initializing
% We need to choose the particle's initial phase space position.

% Initialize the particle in the center of the trap.
% We multiply by d (the grid spacing) because the integrator expects
%   physical units, not indices.
xx1 = d * double(dimensions(1) + 1) / 2.0; % mm
yy1 = d * double(dimensions(2) + 1) / 2.0; % mm
zz1 = d * double(dimensions(3) + 1) / 2.0; % mm

% If you are curious to see what the potential from electrode1 would be, 
% compared to the starting point:
potential_maps_reshaped = reshape(potential_maps, [length(electrode_names) dimensions]);

vxxs = zeros([1 n_particles]);
vyys = zeros([1 n_particles]);
vzzs = zeros([1 n_particles]);
xxs  = xx1 + normrnd(0, cloud_radius, [1 n_particles]);
yys  = yy1 + normrnd(0, cloud_radius, [1 n_particles]);
zzs  = zz1 + normrnd(0, cloud_length, [1 n_particles]);

dissoc_speed = 0.5; %mm / us
for idx = 1:n_particles
    % If we put the ion at the center of the trap with zero initial speed, it
    %   will just sit there quietly.  A nice way to give it some speed is by
    %   specifying a temperature and choosing a speed based on that.
    maxwell = @(v) maxwell_pdf(v, m, T);
    % Maxwell-Boltzmann distribution speeds
    v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
    % Uniform distribution direcitions
    theta = 2 * pi * rand();
    phi = 2 * pi * rand();
    % Convert speed & direction -> velocity
    vxxs(idx) = v * sin(theta) * cos(phi);      % mm / us
    vyys(idx) = v * sin(theta) * sin(phi);      % mm / us
    vzzs(idx) = v * cos(theta);                 % mm / us

    vxxs(idx) = vxxs(idx) + dissoc_speed * (-1)^idx;
end

cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);
vzzs = vzzs - 1.0e-3 * sqrt(2 * axial_eV * cmr) * ones(size(vzzs));

%% Integration: FORTRAN version
tic
potential_maps_size = size(potential_maps);
potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
[xss, yss, zss, tss, itss] = ...
    fly_ensemble(16384, n_particles, xxs, yys, zzs, vxxs, vyys, vzzs, ...
                 potential_maps, voltages, step_times, ...
                 time_steps, dimensions, int32(is_electrode), length(electrode_names), m, q, d, ...
                 maxdist, end_time);
elapsed_time = toc;
fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(sum(itss) / elapsed_time));


% %%
% 
% plot_idx = 3;
% plot3(xss(plot_idx, :), yss(plot_idx, :), zss(plot_idx, :))
% axis image
% 
% 
% % xlim([min(tss(plot_idx, :)) max(tss(plot_idx, :))])

%%

fh = figure(...
    Name            = '',...
    WindowStyle     = 'docked',...
    WindowState     = 'normal',...
    Color           = 'w',...
    Units           = 'pixel',...
    ... Position        = [0 0 args.HeightResolution*args.AspectRatio args.HeightResolution],...
    PaperUnits      = 'centimeters',...
    PaperSize       = [9.0 9.0/((1+sqrt(5))/2)],...
    PaperPosition   = [0 0 9.0 9.0/((1+sqrt(5))/2)],...
    Resize          = 'on'...
    );

%% Projection of electrode onto surface.
% Looking along the x direction

y_center = d * double(dimensions(2) + 1) / 2.0;
z_center = d * double(dimensions(3) + 1) / 2.0;
ys = d*(1:dimensions(2)) - y_center;
zs = d*(1:dimensions(3)) - z_center;
cut = 25;
projected = ~(squeeze(sum(is_electrode(cut:(end-cut), :, :), 1)) > 0);
imagesc(zs, ys, projected)
colormap('gray')
axis image
xlim([-50 50])
ylim([-20 20])

% Plot trajectories
hold on
for idx = 1:n_particles
    plot(zss(idx, :) - z_center, yss(idx, :) - y_center, '-', Color="#1E2F97");
end
hold off


%%

ax = gca;
ax.PlotBoxAspectRatio      = [((1+sqrt(5))/2) 1 1];
ax.Box                     = 'on';
ax.TickLabelInterpreter    = 'latex';
ax.FontName                = 'TimesNewRoman';
ax.FontSize                = 12;
ax.XMinorTick              = 'on';
ax.YMinorTick              = 'off';
ax.XGrid                   = 'off';
ax.YGrid                   = 'on';



