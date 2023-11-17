% Variables to set beforehand

simion_path = "SIM001_data/001_linear_paul_trap";

electrode_names = ["quad_doubled.PA1.patxt", ...
                   "quad_doubled.PA2.patxt", ...
                   "quad_doubled.PA3.patxt", ...
                   "quad_doubled.PA4.patxt", ...
                   "quad_doubled.PA5.patxt", ...
                   "quad_doubled.PA6.patxt", ...
                   "quad_doubled.PA7.patxt", ...
                   "quad_doubled.PA8.patxt", ...
                  ];

start_line = 22;

d = 1.0;

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

potential = reshape(tensorprod(potential_maps, [1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0], 1, 2), dimensions);
% Take a slice at z/2
potential = squeeze(potential(:, :, round(dimensions(3)/2)));
electrode_slice = squeeze(is_electrode(:, :, round(dimensions(3)/2)));

%% Plot
imagesc(potential)
axis image

%% Rejigger the current slice at z/2 into a polar coordinate system.

origin = double(dimensions(1:2)) / 2.0;
r_max = min(origin) - d;

% O is theta
x = @(r, O) r .* cos(O);
y = @(r, O) r .* sin(O);
potential = @(r, O) linInterpolate_many2D(potential, x(r, O) + origin(1), y(r, O) + origin(2), d);

res = 256;
rs = zeros([1, res^2]);
Os = zeros([1, res^2]);
r_spots = linspace(0, r_max, res);
O_spots = linspace(0, 2 * pi, res);
for idx1 = 1:res
    for idx2 = 1:res
        rs((idx1 - 1) * res + idx2) = r_spots(idx1);
        Os((idx1 - 1) * res + idx2) = O_spots(idx2);
    end
end
figure
polarscatter(Os, rs, 3, potential(rs, Os));

%% Find radius for source-free solution

r_min = r_max;
for x = 1:dimensions(1)
    for y = 1:dimensions(2)
        r = sqrt((double(x) - origin(1))^2.0 + (double(y) - origin(2))^2.0);
        if electrode_slice(x,y) && (r < r_min)
            r_min = r;
        end
    end
end

%% Potentials along various directions

figure
rs = linspace(0, r_min, 1024);
hold on
for O = [0, pi/4, pi/2]
    plot(rs, potential(rs, repmat(O, [1 1024])));
end
hold off

%% Functions

function results = linInterpolate_many2D(potential, xs, ys, d)
    results = zeros([1 length(xs)]);
    for idx = 1:length(xs)
        results(idx) = linInterpolate2D(potential, xs(idx), ys(idx), d);
    end
end

