close all
electrode_names = ["comparison_pa.pa1.patxt", ...
                   "comparison_pa.pa2.patxt"
                  ];
start_line = 19;
loadanyways = false;
[potential_maps, is_electrode, dimensions] = ...
    readFile(electrode_names, start_line);
potential_maps = potential_maps / 10000.0;
is_electrode = logical(is_electrode);

map = potential_maps(1, :) - potential_maps(2, :);

map = reshape(map, [1024 1024])';

imagesc(map)
set(gca,'YDir','normal')
drawnow;

%%

figure

start = 4;
final = 10;
final_accuracy = 8.0;
d_final = 0.1682;

d_start = d_final * 2.0^(final - start);


refinement_number = 0;
for res_exp = start:final
    refinement_number = refinement_number + 1;
    res = 2^res_exp;
    dimensions = [res res];
    if refinement_number == 1
        potential_array = zeros(dimensions ./ 2);
    end
    boundary_conditions = zeros(dimensions);

    % Refined boundary conditions
    for r = 1:dimensions(1)
        for z = 1:dimensions(2)
            % if r == 2 * round(res / 10.0)
            %     if z > 2 * round(res / 5.0) && z < 3 * round(res / 5.0)
            %         boundary_conditions(r, z) = 10.0;
            %     end
            % end
            % if r == 4 * round(res / 10.0)
            %     boundary_conditions(r, z) = 30.0;
            % end
            % if r == 8 * round(res / 10.0)
            %     if z > round(res / 4.0) && z < 3 * round(res / 4.0)
            %         boundary_conditions(r, z) = -10.0;
            %     end
            % end
            % if (r - 14.0 * double(dimensions(1)) / 29.0)^2.0 ...
            %         + (z - 14.0 * double(dimensions(1)) / 29.0)^2.0 ...
            %             < double(dimensions(1)) / 32
            %     boundary_conditions(r, z) = 1.0;
            % elseif (r - 15.0 * double(dimensions(1)) / 29.0)^2.0 ...
            %         + (z - 15.0 * double(dimensions(1)) / 29.0)^2.0 ...
            %             < double(dimensions(1)) / 32
            %     boundary_conditions(r, z) = -1.0;
            % end
            if z >= round(double(dimensions(2)) / 8.0) && z <= round(7.0 * double(dimensions(2)) / 8.0)
                if r == round(double(dimensions(1)) / 8.0)
                    boundary_conditions(r, z) =  1.0;
                elseif r == round(2.0 * double(dimensions(1)) / 8.0)
                    boundary_conditions(r, z) = -1.0;
                end
            end
        end
    end

    imagesc(potential_array)
    set(gca,'YDir','normal')
    drawnow;
    
    % Refine potential array
    new_potential_array = zeros(dimensions);
    for r = 2:res + 1
        for z = 2:res + 1
            % Crude 'nearest-neighbor' interpolation.
            new_potential_array(r - 1, z - 1) ...
                = potential_array(floor(r / 2), floor(z / 2));
        end
    end

    potential_array = solve_laplace(min([d_start * double(2.0^(start - res_exp)), double(res / 64.0)]), new_potential_array, boundary_conditions, 8);
end

imagesc(potential_array)
set(gca,'YDir','normal')
drawnow;

% Final refinement
potential_array = solve_laplace(d_final, potential_array, boundary_conditions, final_accuracy);

imagesc(potential_array)
set(gca,'YDir','normal')

figure
imagesc(abs(map - potential_array))
title(sum(abs(map(:) - potential_array(:))) / (1024 * 1024))

%%
function solved = solve_laplace(d, potential_array, boundary_conditions, accuracy)
    dimensions = size(potential_array);
    for r = 2:dimensions(1)
        for z = 2:dimensions(2)
            if boundary_conditions(r, z) ~= 0.0
                potential_array(r, z) = boundary_conditions(r, z);
            end
        end
    end
    
    change = 0.0;
    last_change = 0.0;
    i = 0;
    while (last_change - change) / (dimensions(1) * dimensions(2)) > 10.0^(-accuracy) ...
            || i < 64 %min([i < 1 * max(dimensions), 256])
        i = i + 1;
        last_change = change;
        change = 0.0;
        for r = 2:(dimensions(1) - 1)
            potential_array(r, end) = potential_array(r, end - 1);
            potential_array(r, 1) = potential_array(r, 2);
        end
        for z = 2:(dimensions(2) - 1)
            potential_array(end, z) = potential_array(end - 1, z);
            potential_array(1, z) = potential_array(2, z);
        end
        for z = 2:dimensions(2) - 1
            for r = 2:dimensions(1) - 1
                % Non-boundaries
                if boundary_conditions(r, z) == 0.0
                    prev = potential_array(r, z);
                    potential_array(r, z) = (1 / 4.0) * ( ...
                         (d / (2.0 * r)) * (potential_array(r + 1, z) - potential_array(r - 1, z)) ...
                       + (potential_array(r + 1, z) + potential_array(r - 1, z)) ...
                       + (potential_array(r, z + 1) + potential_array(r, z - 1)));
                    change = change + abs(potential_array(r, z) - prev);
                end
                potential_array(1, z) = (1 / 4.0) * (2 * potential_array(2, z) ...
                       + potential_array(1, z + 1) + potential_array(1, z - 1));
            end
        end
        for z = (dimensions(2) - 1):-1:2
            for r = (dimensions(1) - 1):-1:2
                % Non-boundaries
                if boundary_conditions(r, z) == 0.0
                    prev = potential_array(r, z);
                    potential_array(r, z) = (1 / 4.0) * ( ...
                         (d / (2.0 * r)) * (potential_array(r + 1, z) - potential_array(r - 1, z)) ...
                       + (potential_array(r + 1, z) + potential_array(r - 1, z)) ...
                       + (potential_array(r, z + 1) + potential_array(r, z - 1)));
                    change = change + abs(potential_array(r, z) - prev);
                end
                potential_array(1, z) = (1 / 4.0) * (2 * potential_array(2, z) ...
                       + potential_array(1, z + 1) + potential_array(1, z - 1));
            end
        end
        fprintf("%d (%d x %d): %.3g, delta = %.3g\n", ...
            i, dimensions(1), dimensions(2), ...
            change / (dimensions(1) * dimensions(2)), ...
            (last_change - change) / (dimensions(1) * dimensions(2)))
    end
    solved = potential_array;
end