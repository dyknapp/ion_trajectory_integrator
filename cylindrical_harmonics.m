clear;
p = 0;
n_max = 100;
bessel_zeros = besselzero(p, n_max, 1);

% %% Plot Ha cross-section
% 
% is_to_show = 8;
% % figure('Visible', 'off');
% tiledlayout(is_to_show, 1);
% fig = gcf;
% % set(fig, 'papertype', 'A4');
% % set(fig, 'paperorientation', 'portrait');
% for i = 1:is_to_show
%     r0 = 1.0;
%     L  = 2 * pi;
% 
%     res = 1024;
%     theta = pi / 2;
%     rs = linspace(0, r0, res);
%     zs = linspace(0,  L, res);
%     phis = zeros(res);
% 
%     for row_idx = 1:res
%         phis(row_idx, :) = Ha(p, i, i, i, bessel_zeros, r0, L, rs, theta, zs(row_idx));
%     end
% 
%     nexttile;
%     imagesc(zs, rs, phis');
%     axis image; set(gca,'YDir','normal'); colormap('pink');
%     title(sprintf("i = j = k = %d", i));
% end
% 
% % name = 'cylindrical_harmonics_examples';
% % print(fig, name, '-dpdf', '-fillpage')
% % winopen(sprintf('%s.pdf', name))

%%  Plot example PA

% Names of individual files containing electrodes' potentials
electrode_names = ["cylindrical_harmonics_test.pa1.patxt"];

% Where does the data in the patxt file start?
start_line = 19;

loadanyways = true;
if ~exist('dimensions', 'var') || loadanyways
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end
potential_maps = squeeze(reshape(potential_maps, dimensions));
d = 1.0;

% imagesc(d*(1:dimensions(2)), d*(1:dimensions(1)), potential_maps')
% axis image; set(gca,'YDir','normal'); colormap('pink');

% %% Find the expansion coefficients
% 
% % Parameters:
r0 = double(dimensions(1)) * d;
L  = double(dimensions(2)) * d;
% 
% % Regularizer:
% rs = zeros(dimensions);
% for idx = 0:(dimensions(1) - 1)
%     rs(:, idx + 1) = double(idx) / double(dimensions(1) - 1);
% end
% 
% % Harmonic:
% max_max_coeff = 15;
% figure('visible', 'off');
% tiledlayout(5, 3);
% 
% for max_coeff = 1:max_max_coeff
%     tic
%     coefficients = zeros([max_coeff max_coeff max_coeff]);
%     parfor i = 1:max_coeff
%         for j = 1:max_coeff
%             for k = 1:max_coeff
%                 theta = pi / 2;
%                 rs = linspace(0, r0, dimensions(1));
%                 zs = linspace(0,  L, dimensions(2));
%                 phis = zeros(dimensions);
%                 for row_idx = 1:dimensions(1)
%                     phis(row_idx, :) = Ha(p, i, j, k, bessel_zeros, r0, L, rs, theta, zs(row_idx));
%                 end
% 
%                 % Take the product
%                 product = rs .* phis .* potential_maps;
%                 coefficients(i, j, k) = sum(product, "all") * 2 / (pi * besselj(p + 1, bessel_zeros(i)))^2;
%             end
%         end
%     end
% 
%     % Reconstruct the potential
%     reconstructed_potential = zeros(dimensions);
%     parfor i = 1:max_coeff
%         for j = 1:max_coeff
%             for k = 1:max_coeff
%                 theta = pi / 2;
%                 rs = linspace(0, r0, dimensions(1));
%                 zs = linspace(0,  L, dimensions(2));
%                 phis = zeros(dimensions);
%                 for row_idx = 1:dimensions(1)
%                     phis(row_idx, :) = Ha(p, i, j, k, bessel_zeros, r0, L, rs, theta, zs(row_idx));
%                 end
% 
%                 reconstructed_potential = reconstructed_potential  ...
%                     + coefficients(i, j, k) * phis;
%             end
%         end
%     end
% 
%     % Cheeky normalization
%     reconstructed_potential = reconstructed_potential / max(reconstructed_potential, [], "all");
% 
%     % plot
%     nexttile
%     imagesc(reconstructed_potential')
%     axis image; set(gca,'YDir','normal'); colormap('pink');
%     title(sprintf('i,j,k < %d', max_coeff));
%     toc
%     drawnow('update');
% end
% 
% name = 'cylindrical_harmonics_terms';
% fig = gcf;
% set(fig, 'papertype', 'A4');
% set(fig, 'paperorientation', 'portrait');
% print(fig, name, '-dpdf', '-fillpage')
% winopen(sprintf('%s.pdf', name))
% clf

%%

max_coeff = 32;
coefficients = zeros([max_coeff max_coeff max_coeff]);
parfor i = 1:max_coeff
    for j = 1:max_coeff
        for k = 1:max_coeff
            theta = pi / 2;
            rs = linspace(0, r0, dimensions(1));
            zs = linspace(0,  L, dimensions(2));
            phis = zeros(dimensions);
            for row_idx = 1:dimensions(1)
                phis(row_idx, :) = Ha(p, i, j, k, bessel_zeros, r0, L, rs, theta, zs(row_idx));
            end

            % Take the product
            product = rs .* phis .* potential_maps;
            coefficients(i, j, k) = sum(product, "all") * 2 / (pi * besselj(p + 1, bessel_zeros(i)))^2;
        end
    end
end

% Reconstruct the potential
reconstructed_potential = zeros(dimensions);
parfor i = 1:max_coeff
    for j = 1:max_coeff
        for k = 1:max_coeff
            theta = pi / 2;
            rs = linspace(0, r0, dimensions(1));
            zs = linspace(0,  L, dimensions(2));
            phis = zeros(dimensions);
            for row_idx = 1:dimensions(1)
                phis(row_idx, :) = Ha(p, i, j, k, bessel_zeros, r0, L, rs, theta, zs(row_idx));
            end

            reconstructed_potential = reconstructed_potential  ...
                + coefficients(i, j, k) * phis;
        end
    end
end

% Cheeky normalization
reconstructed_potential = reconstructed_potential / max(reconstructed_potential, [], "all");

% plot
imagesc(reconstructed_potential')
axis image; set(gca,'YDir','normal'); colormap('pink');
title(sprintf('i,j,k < %d', max_coeff));

figure
error_T = (reconstructed_potential' - potential_maps');% ./ potential_maps';
error_T = error_T .* (is_electrode == 0)';
imagesc(error_T)
axis image; set(gca,'YDir','normal'); colormap('pink');
title(sprintf('i,j,k < %d, RMS error = %.2g', max_coeff, sqrt(mean(error_T.^2.0, "all")))); colorbar();
shg

%% Functions

function phi = Ha(p, i, j, k, bessel_zeros, r0, L, r, theta, z)
    arguments
        p            (1, 1) int32
        i            (1, 1) int32
        j            (1, 1) int32
        k            (1, 1) int32
        bessel_zeros (1, :) {double, mustBeLEQ_BZ(bessel_zeros, i)}
        r0           (1, 1) double
        L            (1, 1) double
        r            (1, :) double
        theta        (1, :) double
        z            (1, :) {double, oneArray(r, theta, z)}
    end
    phi = sin(double(k) * pi * z / L);
    phi = phi .* cos(double(j - 1) * theta);
    phi = phi .* besselj(double(p), bessel_zeros(i) * r / r0);
end

function phi = Hb(p, i, j, k, bessel_zeros, r0, L, r, theta, z)
    arguments
        p            (1, 1) int32
        i            (1, 1) int32
        j            (1, 1) int32
        k            (1, 1) int32
        bessel_zeros (1, :) {double, mustBeLEQ_BZ(bessel_zeros, i)}
        r0           (1, 1) double
        L            (1, 1) double
        r            (1, :) double
        theta        (1, :) double
        z            (1, :) {double, oneArray(r, theta, z)}
    end
    phi = sin(double(k) * pi * z / L);
    phi = phi .* cos(double(j - 1) * theta);
    phi = phi .* besselj(double(p), bessel_zeros(i) * r / r0);
end

function mustBeLEQ_BZ(a, b)
    if length(a) < b
        eid = 'Length:TooShort';
        msg = 'There is insufficient data in bessel_zeros for inputted n.';
        error(eid, msg);
    end
end

function oneArray(a,b,c)
    if sum([length(a) length(b) length(c)] == 1) ~= 2
        eid = 'Coords:TooManyVectors';
        msg = 'Input can be vectorized only along one dimension at a time.';
        error(eid, msg);
    end
end