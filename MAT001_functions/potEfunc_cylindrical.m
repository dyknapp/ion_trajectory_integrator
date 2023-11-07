function [E_r, E_z] = potEfunc_cylindrical(r, z, potential, d)
% function to get the electric field from the potential
%
% dknapp, 10.31.2023
%
%
% INPUT:
%       r, z         ... floating point coordinates, not indices, at which matrix needs to be interpolated
%       potential       ... 3D array of potential values in V
%       d               ... physical distance between matrix gridpoints
    

% OUTPUT:
%       E_r, E_z   ... returns E in the x, y, z directions in []
    
    
    E_r = -(linInterpolate2D(potential, abs(r+d/2), z, d) ...
                - linInterpolate2D(potential, abs(r-d/2), z, d)) / (d);
    E_z = -(linInterpolate2D(potential, abs(r), z+d/2, d) ...
                - linInterpolate2D(potential, abs(r), z-d/2, d)) / (d);
end


