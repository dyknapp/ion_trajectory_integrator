function [matrix_interp] = linInterpolate2D(matrix, r, z, d)
% function to interpolate any matrix in 2D linearly
% NOTE: this function is specifically for integrate_trajecoty_cylindrical.
% it adds an offset to the r to accomodate the array index convention.
%
% dknapp, 31.10.2023
%
%
% INPUT:
%       matrix          ... 2D array to interpolate between
%       r, z         ... floating point coordinates, not indices, at which matrix needs to be interpolated
%       d               ... physical distance between matrix gridpoints
    

% OUTPUT:
%       matrix_interp   ... floating point interpolated value   

    % Rescale to match lattice
    r = r / d;
    z = z / d;

    % find the index which represents the coordinates
    r_grid = fix(r);
    z_grid = fix(z);

    % get the coordinate relative to the calculated index
    r_rel = r - r_grid;
    z_rel = z - z_grid;

    % store 8 of the nearest gridpoints to interpolate between
    pot_a = matrix(r_grid, z_grid); 
    pot_b = matrix(r_grid + 1, z_grid);
    pot_d = matrix(r_grid, z_grid + 1);
    pot_f = matrix(r_grid + 1, z_grid + 1); 
  
    % interpolate matrix along the gridlines
    pot_ad = z_rel * pot_d + (1 - z_rel) * pot_a;
    pot_bf = z_rel * pot_f + (1 - z_rel) * pot_b;

    % interpolate matrix along faces of the evaluated grid

    % get the interpolated value
    matrix_interp = r_rel * pot_bf + (1 - r_rel) * pot_ad;
end