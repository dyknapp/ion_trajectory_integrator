% Variables to set beforehand

simion_path = "SIM001_data/001_linear_paul_trap";

electrode_names = ["essentials_only_rf.patxt", ...
                   "essentials_only_dc1.patxt", ...
                   "essentials_only_dc2.patxt", ...
                   "essentials_only_dc3.patxt", ...
                   "essentials_only_dc4.patxt", ...
                   "essentials_only_dc5.patxt", ...
                  ];

start_line = 19;

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

d = 1.0;

m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

% UNIT CONVERSIONS
constants = physical_constants();
m = m * constants("atomic mass unit");
q = q * constants("elementary charge");

%% RF and endcap electrode voltages.

RF_frequency = 40.0e+6;
RF_amplitude = 500.0;
endcap_voltage = 10.0;

%% Plot pseudopotential.

potential_maps = reshape(potential_maps, [length(electrode_names) dimensions]);

pseudo = zeros(dimensions);
% RF pseudopotential
[Ex, Ey, Ez] = gradient(squeeze(potential_maps(1, :, :, :)));
Esquared = (Ex.^2.0 + Ey.^2.0 + Ez.^2.0) / (d / 1.0e+3)^4.0;
pseudo = pseudo + (q^2.0 / (4 * m * (2 * pi * RF_frequency)^2.0)) * Esquared;

% Convert to eV
pseudo = pseudo / constants("elementary charge");
cscale = [0 max(pseudo(:))];

res_mult = 8;
X = (0:double(dimensions(1)-1)) * d;
Y = (0:double(dimensions(2)-1)) * d;
Z = (0:double(dimensions(3)-1)) * d;
Xq = linspace(0, double(dimensions(1)-1) * d, dimensions(1) * res_mult);
Yq = linspace(0, double(dimensions(2)-1) * d, dimensions(2) * res_mult);
Zq = linspace(0, double(dimensions(3)-1) * d, dimensions(3) * res_mult);
[Xin, Yin, Zin] = meshgrid(Xq, Yq, Zq);
pseudo = interp3(X, Y, Z, pseudo, Xin, Yin, Zin);

t = tiledlayout(1, 3);
nexttile
imagesc(Yq, Zq, squeeze(pseudo(round(dimensions(1) * res_mult/2), :, :))');
axis image
colormap turbo
caxis(cscale);

nexttile
imagesc(Xq, Zq, squeeze(pseudo(:, round(dimensions(2) * res_mult/2), :))');
axis image
colormap turbo
caxis(cscale);

nexttile
imagesc(Xq, Yq, squeeze(pseudo(:, :, round(dimensions(3) * res_mult/2)))');
axis image
colormap turbo
caxis(cscale);
colorbar
