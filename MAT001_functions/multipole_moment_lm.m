% Calculate the (l, m) multipole moment
%
% 12.11.2023,  dknapp: wrote the function

function moment = multipole_moment_lm(l, m, x0, y0, z0, d, potential)
    arguments
        % ls (1, :) int32
        % ms (1, :) int32
        l int32
        m int32
        x0 double
        y0 double
        z0 double
        d  double
        potential (:, :, :) double
    end
    dims = size(potential);
    % Triple integral:
    dV = d^3;
    integral = 0;
    for x_grid = 1:dims(1)
        for y_grid = 1:dims(2)
            for z_grid = 1:dims(3)
                x = x_grid - x0;
                y = y_grid - y0;
                z = z_grid - z0;
                
                [phi, theta, r] = cart2sph(x, y, z);
                
                integral = integral ...
                    + r^double(l) ...
                    * potential(x_grid, y_grid, z_grid) ...
                    * Ylm(l, m, theta, phi) ...
                    * dV;
            end
        end
    end

    moment = integral;
end