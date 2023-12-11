constants = physical_constants();

simion_path = "SIM001_data/008_quadrupole_deflector";

% Names of individual files containing electrodes' potentials
electrode_names = ["COLQUAD5_modified.PA1.patxt", ...
                   "COLQUAD5_modified.PA2.patxt", ...
                   "COLQUAD5_modified.PA3.patxt", ...
                   "COLQUAD5_modified.PA4.patxt", ...
                   "COLQUAD5_modified.PA5.patxt", ...
                   "COLQUAD5_modified.PA6.patxt", ...
                   "COLQUAD5_modified.PA7.patxt", ...
                   "COLQUAD5_modified.PA8.patxt"
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
voltages(4) = -10.0;

% Set hyperbolas to voltages
voltages(2) =  21.5;
voltages(3) = -21.5;

% shims should have opposite voltages as their corresponding hyperbolas
voltages(7) =  1.0;
voltages(8) = -1.0;

%%  Set up starting conditions

side_lines = 10;
initial_energy = 10.0;      %eV
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
