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

% Set hyperbolas to voltages
voltages(2) =  23.5;
voltages(3) = -23.5;

% shims should have opposite voltages as their corresponding hyperbolas
voltages(7) =  1.0;
voltages(8) = -1.0;

%%  Set up starting conditions

initial_energy = 10.0;      %eV
initial_angle = pi / 4.0;
initial_location = [5 * d, 5 * d]; %mm

xx = initial_location(1);
yy = initial_location(2);

v = sqrt(2 * initial_energy * constants("elementary charge") / (m * constants("atomic mass unit"))) * 1.0e-3;
vxx = v * cos(initial_angle);
vyy = v * sin(initial_angle);

[x_traj, y_traj, ts, its] = ray_optics_2D(xx, yy, vxx, vyy, potential_maps, voltages, dimensions, length(electrode_names), m, q, d, maxdist, end_time);

its = int32(its);
x_traj = x_traj(1:its);
y_traj = y_traj(1:its);
ts     =     ts(1:its);

check = squeeze(tensorprod(voltages, potential_maps, 1, 1));
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), check)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("turbo");

%%
figure
imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode)
xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
axis image
hold on
plot(x_traj, y_traj, '-g');
hold off
