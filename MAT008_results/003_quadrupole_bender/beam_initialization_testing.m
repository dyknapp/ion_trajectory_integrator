%% Initialize a blank potential
constants = physical_constants();

electrode_names = ["blank"];
dimensions = int32([512 128]);
potential_maps = reshape(zeros(dimensions), [1 dimensions]);
is_electrode = int32(zeros(dimensions));
voltages = zeros(size(electrode_names));

d = 1.0;                 % mm/gu
start_time =    0.0;     % us
end_time   =  100.0;     % us
m = 2.0;                 % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                 % atomic units
maxdist =  0.00025;      % mm
sample_dist = 0.5;       % mm

%%  Set up starting conditions

collimation_f = 150.;
side_lines = 10;
initial_energy = 1.0;      %eV
initial_angle = 0.0;
initial_location = double([3.*d dimensions(2)/2.]); %mm
initial_spread = double(dimensions(2))/2.; %mm

xx = initial_location(1);
yy = initial_location(2);

offsets = double(-side_lines:side_lines) * initial_spread * 0.5 / double(side_lines);
particles = length(offsets);

xxs = zeros([particles 1]);
yys = zeros([particles 1]);

v = sqrt(2 * initial_energy * constants("elementary charge") / (m * constants("atomic mass unit"))) * 1.0e-3;
for idx = 1:particles
    xxs(idx) = xx + offsets(idx) * -sin(initial_angle);
    yys(idx) = yy + offsets(idx) *  cos(initial_angle);
    
    focus_angle = atan2(offsets(idx), collimation_f);

    vxxs(idx) = v * cos(initial_angle - focus_angle);
    vyys(idx) = v * sin(initial_angle - focus_angle);
end

%% Run the trajectories

% trajectories = zeros([length(offsets), 1024, 3]);
% deaths = zeros([1 length(offsets)]);
% itss = deaths; datass = deaths;
% for idx = 1:length(offsets)
%     position = [ xxs(idx)  yys(idx)];
%     velocity = [vxxs(idx) vyys(idx)];
%     [trajectories(idx, :, :), deaths(idx), itss(idx), datass(idx)] = ...
%         ray_optics_spaced(position, velocity, sample_dist, int32(is_electrode), ...
%                           potential_maps, voltages, int32(dimensions), ...
%                           length(electrode_names), m, q, d, maxdist, end_time);
% end

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

hold on
for idx = 1:length(offsets)
    path = squeeze(trajectories(idx, 1:datass(idx), :));
    plot(path(:, 1), path(:, 2), '-k')
end
hold off
axis image
xlim([0 dimensions(1)])
ylim([0 dimensions(2)])

% Where the focus should be
xline(collimation_f + initial_location(1), '--b')





