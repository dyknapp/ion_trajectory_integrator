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

d = 0.284;

start_time =    0.0;      % us
end_time   =  100.0;      % us

m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  0.00025;      % mm

%% Collapse into 2D potentials

dimensions = dimensions(1:2);
potential_maps = reshape(potential_maps, [length(electrode_names) dimensions]);
potential_maps = squeeze(potential_maps);

%% Apply the required mirroring

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

%% View electrode

% % view = 7;
% % 
% % imagesc(squeeze(potential_maps(view, :, :)))
% % axis image

% 1: outer box
% 2: y-direction hyperbolas
% 3: x-direction hyperbolas
% 4: middle einzel element
% 5: outer  einzel element
% 6: middle einzel element, lower
% 7: y-shims ( + lower einzel why?????)  anyways, if we stay upper it
% doesn't matter
% 8: x-shims

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

side_lines = 10;
initial_energy = 1.0;      %eV
initial_angle = pi / 4.0;
initial_location = [15 * d, 15 * d]; %mm
initial_spread = 8; %mm

xx = initial_location(1);
yy = initial_location(2);

offsets = double(-side_lines:side_lines) * initial_spread * 0.5 / double(side_lines);
particles = length(offsets);

xxs = zeros([particles 1]);
yys = zeros([particles 1]);

for idx = 1:particles
    xxs(idx) = xx + offsets(idx) * -sin(initial_angle);
    yys(idx) = yy + offsets(idx) *  cos(initial_angle);
end

v = sqrt(2 * initial_energy * constants("elementary charge") / (m * constants("atomic mass unit"))) * 1.0e-3;
vxxs = v * cos(initial_angle) * ones([particles 1]);
vyys = v * sin(initial_angle) * ones([particles 1]);

ms = m * ones([particles 1]);
qs = q * ones([particles 1]);

voltagess = repmat(voltages, [1, particles])';

[x_trajs, y_trajs, tss, itss, datass] = ray_optics_ensemble(xxs, yys, vxxs, vyys, potential_maps, voltagess, dimensions, length(electrode_names), ms, qs, d, maxdist, end_time);

check = squeeze(tensorprod(voltages, potential_maps, 1, 1));
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), check)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("turbo");
axis image

%% 

figure
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image
hold on
for idx = 1:particles
    if offsets(idx) == 0
        plot(x_trajs(idx, 1:datass(idx)), y_trajs(idx, 1:datass(idx)), '-r');
    else
        plot(x_trajs(idx, 1:datass(idx)), y_trajs(idx, 1:datass(idx)), '-g', 'LineWidth', 0.1);
    end
end
hold off

%% With new integrator

sample_dist = 0.5;
position = [xx yy];
velocity = v * [cos(initial_angle) sin(initial_angle)];

tic
[trajectory, death, its, datas] = ...
    ray_optics_spaced(position, velocity, sample_dist, int32(is_electrode), ...
                      potential_maps, voltages, int32(dimensions), ...
                      length(electrode_names), m, q, d, maxdist, end_time);
toc
trajectory = trajectory(1:datas-1, :);

% plotting
figure
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image
hold on
if death == 1
    plot(trajectory(:, 1), trajectory(:, 2), '-g')
elseif death == 2
    plot(trajectory(:, 1), trajectory(:, 2), '-r')
end
hold off

%% Faster FORTRAN

positions = [xxs yys]'; velocities = [vxxs; vyys];
tic
[trajectories, deaths, itss, datass] = ...
    ray_optics_spaced_ensemble(length(offsets), positions, velocities, ...
                               sample_dist, int32(is_electrode), ...
                               potential_maps, voltages, int32(dimensions), ...
                               length(electrode_names), m, q, d, maxdist, end_time);
toc
trajectories = reshape(trajectories, [length(offsets) 1024 3]);


%% Plot the paths

figure
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image
hold on
for idx = 1:length(offsets)
    if offsets(idx) == 0
        color = 'r';
    else
        color = 'g';
    end
    path = squeeze(trajectories(idx, 1:datass(idx), :));
    plot(path(:, 1), path(:, 2), color)
end
hold off
axis image
xlim([0 dimensions(1)] * d)
ylim([0 dimensions(2)] * d)

%% Liouville's theorem: interpolations and derivatives

% Some stats to help with config
times = squeeze(trajectories(:, :, 3));
fprintf("Time that the last particle dies: %.3gus\n", max(times, [], "all"))
fprintf("Mean death time of all particles: %.3gus\n", mean(max(times, [], 2)))

% For Liouville plot
% 0, 3.5, 8
t = 17;
position_axis = 2;
velocity_axis = 1;
max_shown_displacement = 5.0;
max_shown_speed        = 5.0;

gif = false;
if gif
    !del quadrupole_Liouville.gif
    figure('Position', [0 0 1000 500])
end
for gif_time = linspace(0, mean(max(times, [], 2)), 300)
    if gif
        t = gif_time;
    end
    
    % Construct the positions and velocities
    positions = zeros([1 length(offsets)]); velocities = positions;
    death_times = max(times, [], 2);
    for idx = 1:length(offsets)
        if death_times(idx) >= t
            positions(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [position_axis 3]))',   t);
            velocities(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [velocity_axis 3]))', t);
        end
    end
    % Plot it
    subplot(1, 2, 1); % The phase space
    scatter(positions - mean(positions), velocities - mean(velocities), '.k');
    xlabel(sprintf('position axis %d (mm)', position_axis));
    ylabel(sprintf('velocity axis %d (mm/us)', velocity_axis));
    xlim([-1 1] * max_shown_displacement);
    ylim([-1 1] * max_shown_speed);
    
    subplot(1, 2, 2); % Where are the particles?
    imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
    axis image
    xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
    hold on
    for idx = 1:length(offsets)
        if death_times(idx) >= t
            xs(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [1 3]))',   t);
            ys(idx) = interpolate_position(squeeze(trajectories(idx, 1:datass(idx), ...
                                [2 3]))',   t);
        end
    end
    scatter(xs, ys, '.g', LineWidth=5)
    hold off
    if gif
        exportgraphics(gcf,'quadrupole_Liouville.gif','Append',true);
    else
        break
    end
end

%% FUNCTIONS

function position = interpolate_position(trajectory, time)
    % trajectory should be a (position, time), 2xN row matrix
    position = interp1(trajectory(2, :), trajectory(1, :), time, 'pchip', 'extrap');
end

function velocity = interpolate_velocity(trajectory, time)
    % trajectory should be a (velocity, time), 2xN row matrix
    time_diffs = diff(trajectory(2, :));
    velocities = diff(trajectory(1, :)) / time_diffs;
    middle_times = trajectory(2, 1:(end - 1)) + time_diffs*0.5;
    velocity = interp1(middle_times, velocities, time, 'pchip', 'extrap');
end










