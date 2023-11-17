function [Ex, Ey, Ez] = field_at(x, y, z, d, ...
                                 potential_maps, voltages)
    %% 3D linear interpolation step

    % Find the index which represents the coordinates
    x_grid = floor(x / d);
    y_grid = floor(y / d);
    z_grid = floor(z / d);

    % Cut out the section of the matrix that is necessary to find the field
    potential_maps = potential_maps(:, (x_grid - 1):(x_grid + 2), ...
                                       (y_grid - 1):(y_grid + 2), ...
                                       (z_grid - 1):(z_grid + 2));

    % Since the region has been cut out, we want to similarly offset x,y,z
    % See notes ??? about the indexing here
    % We want the original x,y,z to be mapped to the 2nd of 4 elements
    x = 2 * d + x - x_grid * d;
    y = 2 * d + y - y_grid * d;
    z = 2 * d + z - z_grid * d;
    % By doing this, we are able to just plug the coordinates and potential
    % straight into the existing linInterpolate3D function.

    % Calculate linear sum of electrode contributions.
    potential = tensorprod(voltages, potential_maps, 2, 1);
    potential = reshape(potential, [4, 4, 4]);

    % Find the field based on the potential
    [Ex, Ey, Ez] = potEfunc(x, y, z, potential, d);
end