simion_path = "SIM001_data/001_linear_paul_trap";

% Names of individual files containing electrodes' potentials
electrode_names = ["quad_doubled.PA1.patxt", ...
                   "quad_doubled.PA2.patxt", ...
                  ];

% Where does the data in the patxt file start?
start_line = 22;

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
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

d = 1.0;
potential_maps = reshape(potential_maps, [length(electrode_names) dimensions]);
potential = squeeze(potential_maps(1, :, :, :) + potential_maps(2, :, :, :));

x0 = double(dimensions(1) - 1) / 2.0;
y0 = double(dimensions(2) - 1) / 2.0;
z0 = double(dimensions(3) - 1) / 2.0;

%%

res = 1024;
outer_r = 30.0;
outer_z = 10.0;

%%

thetas = linspace(0, 2*pi, res);
phi_theta = zeros([1 res]);
for idx = 1:res
    phi_theta(idx) = phi(potential, outer_r, thetas(idx), 0.0, x0, y0, z0, d);
end
polarplot(thetas, phi_theta);

%%

phi_theta_z = zeros(res);
parfor idx_t = 1:res
    for idx_z = 1:res
        phi(idx_t, idx_z) = phi(potential, outer_r, thetas(idx), 0.0, x0, y0, z0, d);

%%

function interp = phi(potential, r, theta, z, x0, y0, z0, d)
    x = r * cos(theta) + x0;
    y = r * sin(theta) + y0;
    z = z              + z0;
    interp = linInterpolate3D(potential, x, y, z, d);
end








