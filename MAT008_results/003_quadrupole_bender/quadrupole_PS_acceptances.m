%% Control panel
d = 0.284;                           % mm/gu : 

% Particle setup
% LOTS OF POINTS
particle_position_offsets = 500;
particle_velocity_offsets = 500;
particles = particle_position_offsets * particle_velocity_offsets;

m = 2.0;                             % amu
q = 1.0;                             % electron charges

% Beam setup
initial_energy = 1.0;                % eV
initial_angle =  45. * (pi/180.);    % radians
% collimation_f = 50.;                 % mm      : Distance at which the cloud should be focused
initial_location = [25 * d, 25 * d]; % mm      : Center of cloud

% SET TEMP AND SPREAD HIGH TO GET LOTS OF PHASE SPACE POINTS
initial_spread = 10.0;               % mm      : Standard deviation of Gaussian position initialization
cloud_temperature =  10.0^1;         % Kelvin  : Initial temperature of cloud

% Simulation parameters
end_time   =  100.0;                 % us : Time cutoff for simulation
maxdist =  0.05;                     % mm : Propagation distance for adaptive timestep
sample_dist = 1.0;                   % mm : Interval for recording trajectory data samples

%% Initialize

constants = physical_constants();

simion_path = "SIM001_data/008_quadrupole_deflector";

% Names of individual files containing electrodes' potentials
electrode_names = ["COLQUAD5_modified.PA1.patxt", ...
                   "COLQUAD5_modified.PA3.patxt", ...
                   "COLQUAD5_modified.PA4.patxt", ...
                   "COLQUAD5_modified.PA6.patxt", ...
                   "COLQUAD5_modified.PA7.patxt", ...
                   "COLQUAD5_modified.PA8.patxt", ...
                   "COLQUAD5_modified.PA9.patxt", ...
                   "COLQUAD5_modified.PA10.patxt"
                  ];

start_line = 19;

clear potential_maps is_electrode dimensions
addpath(simion_path)
[potential_maps, is_electrode, dimensions] = ...
    readFile(electrode_names, start_line);
potential_maps = potential_maps / 10000.0;
is_electrode = logical(is_electrode);

%% Collapse into 2D potentials

dimensions = dimensions(1:2);
potential_maps = reshape(potential_maps, [length(electrode_names) dimensions]);
potential_maps = squeeze(potential_maps);

mirror_axis = 2;
temp_dims = dimensions;
temp_dims(mirror_axis) = temp_dims(mirror_axis) * 2 - 1;

temp = zeros([length(electrode_names) temp_dims]);
temp2 = zeros(temp_dims, 'logical');

if mirror_axis == 2
    temp(:, :, (temp_dims(2)/2):end)  = potential_maps;
    temp2(:, (temp_dims(2)/2):end)  = is_electrode;
    temp(:, :, 1:(temp_dims(2)/2))      = flip(potential_maps, 3);
    temp2(:, 1:(temp_dims(2)/2))      = flip(is_electrode, 2);
elseif mirror_axis == 1
    temp(:, (temp_dims(1)/2):end, :)  = potential_maps;
    temp2((temp_dims(1)/2):end, :)  = is_electrode;
    temp(:, 1:(temp_dims(1)/2), :)      = flip(potential_maps, 2);
    temp2(1:(temp_dims(1)/2), :)      = flip(is_electrode, 1);
end

dimensions = temp_dims;
potential_maps = temp;
is_electrode = temp2;

%% Set voltages
voltages = zeros([length(electrode_names) 1]);

% Collimate input
voltages(4) = -1.8;

% Set hyperbolas to voltages
voltages(2) =  2.4;
voltages(3) = -2.4;

% shims should have opposite voltages as their corresponding hyperbolas
voltages(7) =  0.1;
voltages(8) = -0.1;

% Positions for a thermal cloud:
% %%  Set up starting conditions
% 
% % Position spread
% positions = zeros([2 particles]);
% positions(1, :) = normrnd(0.0, initial_spread / sqrt(double(length(initial_location))), [1 particles]);
% positions(2, :) = normrnd(0.0, initial_spread / sqrt(double(length(initial_location))), [1 particles]);
% % Correct for random offset from the center in the positions
% % See at the end for the final offset.  For now, it's more convenient to
% %   stay without including the offset
% positions(1, :) = positions(1, :) - mean(positions(1, :));
% positions(2, :) = positions(2, :) - mean(positions(2, :));
% 
% % The temperature
% velocities = sqrt(2.) * normrnd(0., sqrt(2*constants("boltzmann")*cloud_temperature/(m*constants("proton mass"))), ...
%                             [2, particles]) * 1.0e-3;
% % Correct for random offset in the random velocities
% velocities(1, :) = velocities(1, :) - mean(velocities(1, :));
% velocities(2, :) = velocities(2, :) - mean(velocities(2, :));
% % Add in the bulk movement
% v = sqrt(2*initial_energy*constants("elementary charge")/(m*constants("atomic mass unit"))) * 1.0e-3;
% rotation_matrix = [cos(initial_angle) -sin(initial_angle); ...
%                    sin(initial_angle)  cos(initial_angle)];
% for idx = 1:particles
%     % Calculate the collimation angle:
% 
%     % Rotate distribution to horizontal, temporarily
%     % The inverse of the rotation matrix is the transverse
%     rotated_position = rotation_matrix' * positions(:, idx);
% 
%     % The offset perpendicular to the beam is the y-coordinate
%     % The offset along the beam is the x-coordinate
%     modified_collimation_f = collimation_f - rotated_position(1);
%     collimation_angle = atan2(rotated_position(2), modified_collimation_f);
% 
%     % Step 3: Modify the velocity!
%     velocities(1, idx) = velocities(1, idx) ...
%         + v * cos(initial_angle - collimation_angle);
%     velocities(2, idx) = velocities(2, idx) ...
%         + v * sin(initial_angle - collimation_angle);
% end
% 
% % Add in the position offset
% positions(1, :) = positions(1, :) + initial_location(1);
% positions(2, :) = positions(2, :) + initial_location(2);

%% Initialize a grid in phase space
% Start with particle beam along x.  Then rotate.

max_beam_radius = 10.0;    % mm
max_speed_variation = 3.0; % mm/us

% Calculate initial distributions of positions and velocities
offsets = linspace(-max_beam_radius, max_beam_radius, particle_position_offsets);
v = sqrt(2*initial_energy*constants("elementary charge")/(m*constants("atomic mass unit"))) * 1.0e-3;
velocity_offsets = linspace(v - max_speed_variation, v + max_speed_variation, particle_velocity_offsets);

% Create a position/velocity grid
[p_vals, v_vals] = meshgrid(offsets, velocity_offsets);

positions = zeros([2 particles]);
velocities = zeros([2 particles]);
positions(2, :)  = p_vals(:) + normrnd(0, max_beam_radius / particle_position_offsets, size(p_vals(:)));
velocities(1, :) = v_vals(:) + normrnd(0, max_speed_variation / particle_velocity_offsets, size(v_vals(:)));

positions = [cos(initial_angle) -sin(initial_angle); ...
             sin(initial_angle)  cos(initial_angle)] * positions;
positions(1,:) = positions(1,:) + initial_location(1);
positions(2,:) = positions(2,:) + initial_location(2);
velocities = [cos(initial_angle) -sin(initial_angle); ...
              sin(initial_angle)  cos(initial_angle)] * velocities;

%% Parallel integrate trajectories
parallel_chunks = 4;

tic
idx_borders = round(linspace(1, particles, parallel_chunks + 1));
results = cell([1 parallel_chunks]);
parfor idx = 1:parallel_chunks
    [ts, des, is, das] = ...
        ray_optics_spaced_ensemble_parallel(idx_borders(idx + 1) - idx_borders(idx), ...
                                   positions(:, idx_borders(idx):(idx_borders(idx + 1) - 1)), ...
                                   velocities(:, idx_borders(idx):(idx_borders(idx + 1) - 1)), ...
                                   sample_dist, int32(is_electrode), ...
                                   potential_maps, voltages, int32(dimensions), ...
                                   length(electrode_names), m, q, d, maxdist, end_time);
    results{idx}.trajectories = reshape(ts, [idx_borders(idx + 1) - idx_borders(idx) 1024 3]);
    results{idx}.deaths       = des;
    results{idx}.itss         = is;
    results{idx}.datass       = das;
end

trajectories = zeros([particles 1024 3]);
deaths       = zeros([1 particles]);
itss         = zeros([1 particles]);
datass       = zeros([1 particles]);
for idx = 1:parallel_chunks
    trajectories((idx_borders(idx):(idx_borders(idx + 1) - 1)), :, :) = results{idx}.trajectories;
    deaths((idx_borders(idx):(idx_borders(idx + 1) - 1))) = results{idx}.deaths;
    itss((idx_borders(idx):(idx_borders(idx + 1) - 1))) = results{idx}.itss;
    datass((idx_borders(idx):(idx_borders(idx + 1) - 1))) = results{idx}.datass;
end

toc

%% Plot the paths

% figure
% imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
% xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
% axis image
% hold on
% for idx = 1:particles
%     if datass(idx) > 1
%         if deaths(idx) == 1
%             color = 'g';
%         elseif deaths(idx) == 2
%             color = 'r';
%         else
%             color = 'y';
%         end
%         path = squeeze(trajectories(idx, 1:datass(idx), :));
%         plot(path(:, 1), path(:, 2), color)
%     end
% end
% hold off
% axis image
% xlim([0 dimensions(1)] * d)
% ylim([0 dimensions(2)] * d)

%% Filter out the accepted / rejected coordinates

% Set to > 0.0, just so that we aren't fully trusting pchip extrapolation
t = 0.1;

% Check acceptances
accepted_coords = [];
rejected_coords = [];
for idx = 1:particles
    if datass(idx) >= 4 % We need to have recorded at least a few points for the particle
        if deaths(idx) == 1  ...
                && trajectories(idx, datass(idx), 1) > 100. ...
                    && trajectories(idx, datass(idx), 2) < 40.
    
                accepted_coords(size(accepted_coords, 1) + 1, :) ...
                    = [ interpolate_position(squeeze(trajectories(idx, 1:datass(idx), [1 3]))', t) ...
                        interpolate_position(squeeze(trajectories(idx, 1:datass(idx), [2 3]))', t) ...
                        interpolate_velocity(squeeze(trajectories(idx, 1:datass(idx), [1 3]))', t) ...
                        interpolate_velocity(squeeze(trajectories(idx, 1:datass(idx), [2 3]))', t) ...
                      ];
        else
                rejected_coords(size(rejected_coords, 1) + 1, :) ...
                    = [ interpolate_position(squeeze(trajectories(idx, 1:datass(idx), [1 3]))', t) ...
                        interpolate_position(squeeze(trajectories(idx, 1:datass(idx), [2 3]))', t) ...
                        interpolate_velocity(squeeze(trajectories(idx, 1:datass(idx), [1 3]))', t) ...
                        interpolate_velocity(squeeze(trajectories(idx, 1:datass(idx), [2 3]))', t) ...
                      ];
        end
    end
end

%% Plot them, using the defined beam axis!
% Remember that the x and y we start with are not necessarily useful for
%   thinking about the beam

% Create the linear combinations
accepted_coords(:, 1) = accepted_coords(:, 1) - initial_location(1);
accepted_coords(:, 2) = accepted_coords(:, 2) - initial_location(2);
rejected_coords(:, 1) = rejected_coords(:, 1) - initial_location(1);
rejected_coords(:, 2) = rejected_coords(:, 2) - initial_location(2);

accepted_coords = accepted_coords';
rejected_coords = rejected_coords';

accepted_positions =  [cos(-initial_angle) -sin(-initial_angle); ...
                       sin(-initial_angle)  cos(-initial_angle)] * accepted_coords(1:2, :);
rejected_positions =  [cos(-initial_angle) -sin(-initial_angle); ...
                       sin(-initial_angle)  cos(-initial_angle)] * rejected_coords(1:2, :);
accepted_velocities = [cos(-initial_angle) -sin(-initial_angle); ...
                       sin(-initial_angle)  cos(-initial_angle)] * accepted_coords(3:4, :);
rejected_velocities = [cos(-initial_angle) -sin(-initial_angle); ...
                       sin(-initial_angle)  cos(-initial_angle)] * rejected_coords(3:4, :);

% Perpendicular to input beam
accepted_positions  = accepted_positions(2, :);
rejected_positions  = rejected_positions(2, :);
% Parallel to input beam
accepted_velocities = accepted_velocities(1, :);
rejected_velocities = rejected_velocities(1, :);

%%

figure
hold on
% scatter(accepted_positions - mean(rejected_positions), accepted_velocities - mean(rejected_velocities), '.g');
% scatter(rejected_positions - mean(rejected_positions), rejected_velocities - mean(rejected_velocities), '.r');
scatter(rejected_positions, rejected_velocities, '.r', LineWidth = 0.01);
scatter(accepted_positions, accepted_velocities, '.g', LineWidth = 0.01);
hold off
xlabel('Position perpendicular to beam (mm)')
ylabel('Velocity parallel to beam (mm/us)')


%%
res = 128;

positions_bins = linspace(min(min(accepted_positions), min(rejected_positions)), ...
                          max(max(accepted_positions), max(rejected_positions)), res + 1);
velocities_bins = linspace(min(min(accepted_velocities), min(rejected_velocities)), ...
                          max(max(accepted_velocities), max(rejected_velocities)), res + 1);

h1 = histogram2(rejected_positions, rejected_velocities, ...
    'XBinEdges', positions_bins, 'YBinEdges', velocities_bins);
a_dist = h1.Values;
h2 = histogram2(accepted_positions, accepted_velocities, ...
    'XBinEdges', h1.XBinEdges, 'YBinEdges', h1.YBinEdges);
r_dist = h2.Values;
x_vals = h2.XBinEdges + h2.BinWidth(1);
y_vals = h2.YBinEdges + h2.BinWidth(2);
imagesc(x_vals, y_vals, (a_dist - r_dist)');
% contourf(-a_dist + r_dist, 20)
set(gca,'YDir','normal'); colormap("pink");
xlabel('Position perpendicular to beam (mm)')
ylabel('Velocity parallel to beam (mm/us)')

%% FUNCTIONS
function position = interpolate_position(trajectory, time)
    % trajectory should be a (position, time), 2xN row matrix
    position = interp1(trajectory(2, :), trajectory(1, :), time, 'pchip');
end

function velocity = interpolate_velocity(trajectory, time)
    % trajectory should be a (velocity, time), 2xN row matrix
    time_diffs = diff(trajectory(2, :));
    velocities = diff(trajectory(1, :)) ./ time_diffs;
    middle_times = trajectory(2, 1:(end - 1)) + time_diffs*0.5;
    velocity = interp1(middle_times, velocities, time, 'pchip');
end