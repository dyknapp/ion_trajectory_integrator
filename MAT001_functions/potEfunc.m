function [E_x, E_y, E_z] = potEfunc(x, y, z, potential, d)
% function to get the electric field from the potential
%
% mbroerse, 16.6.2024
%
%
% INPUT:
%       x, y, z         ... floating point coordinates, not indices, at which matrix needs to be interpolated
%       potential       ... 3D array of potential values in V
%       d               ... physical distance between matrix gridpoints
    

% OUTPUT:
%       E_x, E_y, E_z   ... returns E in the x, y, z directions in []
    
    
    E_x = -(linInterpolate3D(potential, x+d/2, y, z, d) ...
                - linInterpolate3D(potential, x-d/2, y, z, d)) / (d);
    E_y = -(linInterpolate3D(potential, x, y+d/2, z, d) ...
                - linInterpolate3D(potential, x, y-d/2, z, d)) / (d);
    E_z = -(linInterpolate3D(potential, x, y, z+d/2, d) ...
                - linInterpolate3D(potential, x, y, z-d/2, d)) / (d);
end


