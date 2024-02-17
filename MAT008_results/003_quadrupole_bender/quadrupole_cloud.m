%% Control panel
d = 0.284;                           % mm/gu : 

% Particle setup
particles = 1000;
m = 2.0;                             % amu
q = 1.0;                             % electron charges

% Beam setup
initial_energy = 1.0;                % eV
initial_angle =  45. * (pi/180.);    % radians
collimation_f = 50.;                 % mm      : Distance at which the cloud should be focused
initial_location = [15 * d, 15 * d]; % mm      : Center of cloud
initial_spread = 5.0;                % mm      : Standard deviation of Gaussian position initialization
cloud_temperature =  10.0;           % Kelvin  : Initial temperature of cloud

% Simulation parameters
end_time   =  100.0;                 % us : Time cutoff for simulation
maxdist =  0.0001;                   % mm : Propagation distance for adaptive timestep
sample_dist = 0.5;                   % mm : Interval for recording trajectory data samples

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
voltages(4) = -1.0;

% Set hyperbolas to voltages
voltages(2) =  2.15;
voltages(3) = -2.15;

% shims should have opposite voltages as their corresponding hyperbolas
voltages(7) =  0.1;
voltages(8) = -0.1;

%%  Set up starting conditions

% Position spread
positions = zeros([2 particles]);
positions(1, :) = normrnd(0.0, initial_spread / sqrt(double(length(initial_location))), [1 particles]);
positions(2, :) = normrnd(0.0, initial_spread / sqrt(double(length(initial_location))), [1 particles]);
% Correct for random offset from the center in the positions
% See at the end for the final offset.  For now, it's more convenient to
%   stay without including the offset
positions(1, :) = positions(1, :) - mean(positions(1, :));
positions(2, :) = positions(2, :) - mean(positions(2, :));

% The temperature
velocities = sqrt(2.) * normrnd(0., sqrt(2*constants("boltzmann")*cloud_temperature/(m*constants("proton mass"))), ...
                            [2, particles]) * 1.0e-3;
% Correct for random offset in the random velocities
velocities(1, :) = velocities(1, :) - mean(velocities(1, :));
velocities(2, :) = velocities(2, :) - mean(velocities(2, :));
% Add in the bulk movement
v = sqrt(2*initial_energy*constants("elementary charge")/(m*constants("atomic mass unit"))) * 1.0e-3;
rotation_matrix = [cos(initial_angle) -sin(initial_angle); ...
                   sin(initial_angle)  cos(initial_angle)];
for idx = 1:particles
    % Calculate the collimation angle:

    % Rotate distribution to horizontal, temporarily
    % The inverse of the rotation matrix is the transverse
    rotated_position = rotation_matrix' * positions(:, idx);

    % The offset perpendicular to the beam is the y-coordinate
    % The offset along the beam is the x-coordinate
    modified_collimation_f = collimation_f - rotated_position(1);
    collimation_angle = atan2(rotated_position(2), modified_collimation_f);

    % Step 3: Modify the velocity!
    velocities(1, idx) = velocities(1, idx) ...
        + v * cos(initial_angle - collimation_angle);
    velocities(2, idx) = velocities(2, idx) ...
        + v * sin(initial_angle - collimation_angle);
end

% Add in the position offset
positions(1, :) = positions(1, :) + initial_location(1);
positions(2, :) = positions(2, :) + initial_location(2);

%% Integrate trajectories

tic
[trajectories, deaths, itss, datass] = ...
    ray_optics_spaced_ensemble(particles, positions, velocities, ...
                               sample_dist, int32(is_electrode), ...
                               potential_maps, voltages, int32(dimensions), ...
                               length(electrode_names), m, q, d, maxdist, end_time);
toc
trajectories = reshape(trajectories, [particles 1024 3]);

% %% Plot the paths
% 
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

%% Liouville's theorem: interpolations and derivatives

% Some stats to help with config
times = squeeze(trajectories(:, :, 3));
fprintf("Time that the last particle dies: %.3gus\n", max(times, [], "all"))
fprintf("Mean death time of all particles: %.3gus\n", mean(max(times, [], 2)))

% For Liouville plot
% 0, 3.5, 8
t = 8;
max_shown_displacement = 10.0;
max_shown_speed        = 10.0;

dims = 2;
subplots_horizontal = 2;
subplots_vertical = 3;

death_times = max(times, [], 2);

gif = false;
if gif
    !del quadrupole_Liouville.gif
    figure('Position', [0 0 750 1000])
end
frame_times = linspace(0, mean(max(times, [], 2)), 300);
emittances = zeros(size(frame_times));
alive_particles = zeros(size(frame_times));
frame_number = 0;
for gif_time = frame_times
    frame_number = frame_number + 1;
    if gif
        t = gif_time;
    end
    
    % Construct the positions and velocities
    positions = zeros([dims particles]);
    velocities = zeros([dims particles]);
    for position_axis = 1:dims
        for velocity_axis = 1:dims
            for idx = 1:particles
                if (death_times(idx) >= t) && (datass(idx) >= 4)
                    positions(position_axis, idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                        [position_axis 3]))',   t);
                    velocities(velocity_axis, idx) = interpolate_velocity(squeeze(trajectories(idx, 1:datass(idx), ...
                                        [velocity_axis 3]))', t);
                else
                    positions(position_axis, idx) = NaN;
                    velocities(velocity_axis, idx) = NaN;
                end
            end
        end
    end
    % Plot it
    subplot_count = 0;
    for position_axis = 1:dims
        for velocity_axis = 1:dims
            subplot_count = subplot_count + 1;
            subplot(subplots_vertical, subplots_horizontal, subplot_count); % The phase space
            idcs = ~(isnan(positions(position_axis, :)) | isnan(velocities(velocity_axis, :)));
            scatter(positions(position_axis, idcs) - mean(positions(position_axis, idcs)), ...
                velocities(velocity_axis, idcs) - mean(velocities(velocity_axis, idcs)), '.k');
            xlabel(sprintf('position axis %d (mm)', position_axis));
            ylabel(sprintf('velocity axis %d (mm/us)', velocity_axis));
            xlim([-1 1] * max_shown_displacement);
            ylim([-1 1] * max_shown_speed);
        end
    end
    
    subplot(subplots_vertical, subplots_horizontal, 5); % Where are the particles?
    imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
    axis image
    xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
    hold on
    for idx = 1:particles
        if (death_times(idx) >= t) && (datass(idx) >= 4)
            xs(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [1 3]))',   t);
            ys(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [2 3]))',   t);
        else
            xs(idx) = NaN;
            ys(idx) = NaN;
        end
    end
    xs = rmmissing(xs);
    ys = rmmissing(ys);
    alive_particles(frame_number) = length(xs);
    scatter(xs, ys, '.g', LineWidth=0.5)
    hold off

    subplot(subplots_vertical, subplots_horizontal, 6);
    imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
    xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
    axis image
    hold on
    for idx = 1:particles
        if datass(idx) > 1
            if deaths(idx) == 1
                color = 'g';
            elseif deaths(idx) == 2
                color = 'r';
            else
                color = 'y';
            end
            path = squeeze(trajectories(idx, 1:datass(idx), :));
            plot(path(:, 1), path(:, 2), color)
        end
    end
    hold off
    axis image
    xlim([0 dimensions(1)] * d)
    ylim([0 dimensions(2)] * d)

    % Phase space volume calculation
    volume = 2 * pi^(dims/2.) / (dims * gamma(dims/2.)) ...
        * sqrt(det(cov([positions' velocities'], "omitrows")));
    emittances(frame_number) = volume;

    sgtitle(sprintf('t = %.3f (us), Phase space volume = %.3g', t, volume))

    if gif
        exportgraphics(gcf,'quadrupole_Liouville.gif','Append',true);
    else
        break
    end
end

%% Emittance variation
figure
yyaxis left
plot(frame_times, emittances);
yyaxis right
plot(frame_times, alive_particles);

xlim([min(frame_times) max(frame_times)])

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


