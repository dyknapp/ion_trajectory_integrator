function [Er, Ez] = field_at_cylindrical(r, z, d, ...
                                 potential_maps, voltages)
    %% 3D linear interpolation step

    % Find the index which represents the coordinates
    % Remember, the r-grid needs to be shifted by 2 to account for 
    % r = -2, -1
    r_grid = (floor(abs(r) / d) + 1) + 2;
    z_grid = (floor(    z  / d) + 1);

    % Cut out the section of the matrix that is necessary to find the field
    potential_maps = potential_maps(:, (r_grid - 2):(r_grid + 2), ...
                                       (z_grid - 2):(z_grid + 2));

    % Since the region has been cut out, we want to similarly offset x,y,z
    % See notes ??? about the indexing here
    % We want the original x,y,z to be mapped to the 3nd of 5 elements
    r = (2 + (r / d) - (r_grid - 1)) * d;
    z = (2 + (z / d) - (z_grid - 1)) * d;
    % By doing this, we are able to just plug the coordinates and potential
    % straight into the existing linInterpolate3D function.

    % Calculate linear sum of electrode contributions.
    potential = tensorprod(voltages, potential_maps, 2, 1);
    potential = reshape(potential, [5, 5]);

    % Find the field based on the potential
    [Er, Ez] = potEfunc_cylindrical(r, z, potential, d);
end