%% Control Panel
constants = physical_constants();

d = 1.0;                           % mm/gu : 

% Particle setup
particles = 1000;
m = 2.0;                             % amu
q = 1.0;                             % electron charges

% Simulation parameters
end_time   =  100.0;                 % us : Time cutoff for simulation
maxdist =  0.001;                    % mm : Propagation distance for adaptive timestep
sample_dist = d;                     % mm : Interval for recording trajectory data samples

%% Initialize

simion_path = "SIM001_data/009_ashfold_stack";

% Names of individual files containing electrodes' potentials
electrode_names = ["ashfold_cross_section.PA2.patxt", ...
                   "ashfold_cross_section.PA4.patxt", ...
                   "ashfold_cross_section.PA5.patxt", ...
                   "ashfold_cross_section.PA6.patxt", ...
                   "ashfold_cross_section.PA10.patxt", ...
                   "ashfold_cross_section.PA11.patxt", ...
                   "ashfold_cross_section.PA13.patxt", ...
                   "ashfold_cross_section.PA14.patxt", ...
                  ];
start_line = 19;

clear potential_maps is_electrode dimensions
addpath(simion_path)
[potential_maps, is_electrode, dimensions] = ...
    readFile(electrode_names, start_line);
potential_maps = potential_maps / 10000.0;
is_electrode = logical(is_electrode);
potential_maps = reshape(potential_maps, [length(electrode_names) dimensions(1) dimensions(2)]);
potential_maps = rot90(potential_maps, [3 2]);

%% View electrodes
electrode_to_view = 3;
potential_maps = squeeze(reshape(potential_maps, [length(electrode_names) dimensions]));
imagesc(squeeze(potential_maps(electrode_to_view, :, :)));
xlabel("z (mm)"); ylabel("r (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image

%% Set voltages:
% Format: 
% [electrode #, time (us), voltage (V);
%  ...        , ...      , ...        ;]
%
% If time is set to zero, the initial condition is defined.
% If there is not specification, the voltage is assumed to be zero.
% When you set a voltage, it stays until the next voltage is set.

% The voltages are set WHEN THE NEXT TIMESTEP WOULD TURN THEM ON.
% The minimum timestep (regardless of max_dist) is end_time * 1.0e-6;
% So the voltage can be turned on as early as t - end_time*1.0e-6

voltages = ...
    [3, 0.1, 1000;];

voltage_lines = size(voltages, 1);
% Make sure that voltages are sorted by time before they get sent to
%   the FORTRAN integration code.

%% Particle Initialization

position = [0, 0, 100];
velocity = [0, 0, 0];


%% Integration

[trajectory, death, its, datas] = ...
    cylindrical_3D_spaced(position, velocity, sample_dist, int32(zeros(size(is_electrode))), ...
                          potential_maps, voltages, voltage_lines, ...
                          dimensions(1:2), length(electrode_names), m, q, d, maxdist, end_time);


%% Plot

plot(trajectory(1:datas, 4), trajectory(1:datas, 3))