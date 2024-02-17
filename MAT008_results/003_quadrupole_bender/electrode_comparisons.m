simion_path = "SIM001_data/008_quadrupole_deflector";

% Names of individual files containing electrodes' potentials
electrode_names = ["COLQUAD5_noshim.PA1.patxt", ...
                   "COLQUAD5_noshim.PA3.patxt", ...
                   "COLQUAD5_noshim.PA4.patxt", ...
                   "COLQUAD5_noshim.PA6.patxt", ...
                   "COLQUAD5_noshim.PA7.patxt", ...
                   "COLQUAD5_noshim.PA8.patxt", ...
                   "COLQUAD5_noshim.PA9.patxt", ...
                   "COLQUAD5_noshim.PA10.patxt", ...
                  ];

start_line = 19;

clear potential_maps is_electrode dimensions
addpath(simion_path)
[potential_maps, is_electrode, dimensions] = ...
    readFile(electrode_names, start_line);
potential_maps = potential_maps / 10000.0;
is_electrode = logical(is_electrode);

d = 0.284;

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

%%

% contourf(squeeze(potential_maps(2, :, :))' - squeeze(potential_maps(3, :, :))', 50)
% axis image


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

%% Integrate trajectories: voltage sweeps
% Set voltages

tests = 20000;
h_voltages = linspace(-5, 5, tests);

acceptances = zeros([1 tests]);
for test_idx = 1:tests
    fprintf("Test %d / %d  :\n", test_idx, tests);

    % Set the voltages for the test
    voltages(2) =  h_voltages(test_idx);
    voltages(3) = -h_voltages(test_idx);

    tic
    [trajectories, deaths, itss, datass] = ...
        ray_optics_spaced_ensemble(particles, positions, velocities, ...
                                   sample_dist, int32(is_electrode), ...
                                   potential_maps, voltages, int32(dimensions), ...
                                   length(electrode_names), m, q, d, maxdist, end_time);
    toc
    trajectories = reshape(trajectories, [particles 1024 3]);
    voltages = zeros([length(electrode_names) 1]);

    % Check acceptances
    for idx = 1:particles
        if deaths(idx) == 1 && datass(idx) >= 1
            % Did we make it out of the right quadrant?
            if trajectories(idx, datass(idx), 1) > 100. ...
                && trajectories(idx, datass(idx), 2) < 40.
                acceptances(test_idx) = acceptances(test_idx) + 1;
            end
        end
    end
end

acceptances = acceptances / particles;

plot(h_voltages, acceptances);
xlabel("Hyperbolic electrode voltage (V)")
ylabel("Acceptance ratio")

%% Plot the paths

[~, test_idx] = max(acceptances);
test_idx = 13202;

% Set the voltages for the test
% Set the voltages for the test
voltages(2) =  h_voltages(test_idx);
voltages(3) = -h_voltages(test_idx);

tic
[trajectories, deaths, itss, datass] = ...
    ray_optics_spaced_ensemble(particles, positions, velocities, ...
                               sample_dist, int32(is_electrode), ...
                               potential_maps, voltages, int32(dimensions), ...
                               length(electrode_names), m, q, d, maxdist, end_time);
toc
trajectories = reshape(trajectories, [particles 1024 3]);

figure
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


