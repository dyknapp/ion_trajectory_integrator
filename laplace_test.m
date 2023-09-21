dimensions = [100 100 100];
potential_array = zeros(dimensions);
boundary_conditions = zeros(dimensions);

d = 1.0;

% Boundary conditions
for x = 1:dimensions(1)
    for y = 1:dimensions(2)
        if   (x - double(dimensions(1) / 2.0))^2 ...
           + (y - double(dimensions(2) / 2.0))^2 ...
           < (double(dimensions(1)) / 5.0)^2
            boundary_conditions(x, y, round(0.4 * dimensions(3))) =  1.0;
            boundary_conditions(x, y, round(0.6 * dimensions(3))) =  1.0;
        end
    end
end

change = dimensions(1) * dimensions(2) * dimensions(3);
i = 0;
while change / (dimensions(1) * dimensions(2) * dimensions(3)) > 2.5e-4
    change = 0;
    i = i + 1;
    % Apply BCs
    for x = 1:dimensions(1)
        for y = 1:dimensions(2)
            for z = 1:dimensions(3)
                if boundary_conditions(x, y, z) ~= 0
                    potential_array(x, y, z) = boundary_conditions(x, y, z);
                end
            end
        end
    end

    % Iteration
    % Apparently the correct operator (DERIVE) is
    for x = 2:dimensions(1) - 1
        for y = 2:dimensions(2) - 1
            for z = 2:dimensions(3) - 1
                prev = potential_array(x, y, z);
                potential_array(x, y, z) = (...
                    + potential_array(x + 1, y    , z    ) ...
                    + potential_array(x - 1, y    , z    ) ...
                    + potential_array(x    , y + 1, z    ) ...
                    + potential_array(x    , y - 1, z    ) ...
                    + potential_array(x    , y    , z + 1) ...
                    + potential_array(x    , y    , z - 1) ...
                    - (d^2) * (0)) / 6;  % Later add charge density
                change = change + abs(potential_array(x, y, z) - prev);
            end
        end
    end
    fprintf("%d: %.3g\n", i, change / (dimensions(1) * dimensions(2) * dimensions(3)))
end
imagesc(squeeze(potential_array(round(dimensions(1) / 2), :, :)))


% % % Compute Laplacian over all space
% % for x = 1:dimensions(1)
% %     for y = 1:dimensions(2)
% %         for z = 1:dimensions(3)
% %             % Laplacian operator
% %             laplace = 0.0;
% %             laplace = laplace + ...
% %                 second_derivative(potential_array(x - 1, y, z), ...
% %                                   potential_array(x    , y, z), ...
% %                                   potential_array(x + 1, y, z), d);
% %             laplace = laplace + ...
% %                 second_derivative(potential_array(x, y - 1, z), ...
% %                                   potential_array(x, y    , z), ...
% %                                   potential_array(x, y + 1, z), d);
% %             laplace = laplace + ...
% %                 second_derivative(potential_array(x, y, z - 1), ...
% %                                   potential_array(x, y, z    ), ...
% %                                   potential_array(x, y, z + 1), d);
% %         end
% %     end
% % end
% % 
% % function dfdt2 = second_derivative(fm, f0, fp, d)
% %     dfdt2 = (fp + fm - 2*f0) / d^2;
% % end