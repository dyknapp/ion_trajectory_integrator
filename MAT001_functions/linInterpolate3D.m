function [matrix_interp] = linInterpolate3D(matrix, x, y, z, d)
% function to interpolate any matrix in 3D linearly
%
% mbroerse, 16.6.2023
% dknapp, 17.8.2023 - reduced number of divisions for speed.
%
%
% INPUT:
%       matrix          ... 3D array to interpolate between
%       x, y, z         ... floating point coordinates, not indices, at which matrix needs to be interpolated
%       d               ... physical distance between matrix gridpoints
    

% OUTPUT:
%       matrix_interp   ... floating point interpolated value   

    % Rescale to match lattice
    x = x / d;
    y = y / d;
    z = z / d;

    % find the index which represents the coordinates
    x_grid = fix(x);
    y_grid = fix(y);
    z_grid = fix(z);

    % get the coordinate relative to the calculated index
    x_rel = x - x_grid;
    y_rel = y - y_grid;
    z_rel = z - z_grid;

    % x_grid = x_grid + 1;
    % y_grid = y_grid + 1;
    % z_grid = z_grid + 1;

    % store 8 of the nearest gridpoints to interpolate between
    pot_a = matrix(x_grid, y_grid, z_grid); 
    pot_b = matrix(x_grid + 1, y_grid, z_grid); 
    pot_c = matrix(x_grid, y_grid + 1, z_grid); 
    pot_d = matrix(x_grid, y_grid, z_grid + 1); 
    pot_e = matrix(x_grid + 1, y_grid + 1, z_grid); 
    pot_f = matrix(x_grid + 1, y_grid, z_grid + 1); 
    pot_g = matrix(x_grid, y_grid + 1, z_grid + 1); 
    pot_h = matrix(x_grid + 1, y_grid + 1, z_grid + 1); 
  
    % interpolate matrix along the gridlines
    pot_cg = z_rel * pot_g + (1 - z_rel) * pot_c;
    pot_ad = z_rel * pot_d + (1 - z_rel) * pot_a;
    pot_eh = z_rel * pot_h + (1 - z_rel) * pot_e;
    pot_bf = z_rel * pot_f + (1 - z_rel) * pot_b;

    % interpolate matrix along faces of the evaluated grid
    pot_x = y_rel * pot_cg + (1 - y_rel) * pot_ad; % a
    pot_xp = y_rel * pot_eh + (1 - y_rel) * pot_bf; % b

    % get the interpolated value
    matrix_interp = x_rel * pot_xp + (1 - x_rel) * pot_x;
end