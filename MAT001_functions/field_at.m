function [Ex, Ey, Ez] = field_at(x, y, z, d, ...
                                 potential_maps, voltages)
    %% 3D linear interpolation step

    % Find the index which represents the coordinates
    x_grid = floor(x / d) + 1;
    y_grid = floor(y / d) + 1;
    z_grid = floor(z / d) + 1;

    % Cut out the section of the matrix that is necessary to find the field
    potential_maps = potential_maps(:, (x_grid - 2):(x_grid + 2), ...
                                       (y_grid - 2):(y_grid + 2), ...
                                       (z_grid - 2):(z_grid + 2));

    % Since the region has been cut out, we want to similarly offset x,y,z
    % See notes ??? about the indexing here
    % We want the original x,y,z to be mapped to the 3nd of 5 elements
    x = (2 + (x / d) - (x_grid - 1)) * d;
    y = (2 + (y / d) - (y_grid - 1)) * d;
    z = (2 + (z / d) - (z_grid - 1)) * d;
    % By doing this, we are able to just plug the coordinates and potential
    % straight into the existing linInterpolate3D function.

    % Calculate linear sum of electrode contributions.
    potential = tensorprod(voltages, potential_maps, 2, 1);
    potential = reshape(potential, [5, 5, 5]);

    % Find the field based on the potential
    [Ex, Ey, Ez] = potEfunc(x, y, z, potential, d);
end