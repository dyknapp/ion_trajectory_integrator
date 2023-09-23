vs = [-15, 100, -5, 100];
rs = [100, 100, 100, 100] ./ 512.0;
z1s = [1, 101, 201, 301] ./ 512.0;
z2s = [90, 190, 290, 400] ./ 512.0;

test(vs, rs, z1s, z2s)

%%

function score = test(vs, rs, z1s, z2s)
    start = 4;
    final = 7;
    final_accuracy = 7.0;
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
    
        % Refined boundary conditions
        [where_bounds, boundary_conditions] = generate_electrodes(vs, rs, z1s, z2s, dimensions);
    
        % imagesc(potential_array)
        % set(gca,'YDir','normal')
        % drawnow;
        
        % Refine potential array
        new_potential_array = zeros(dimensions);
        for r = 2:res + 1
            for z = 2:res + 1
                % Crude 'nearest-neighbor' interpolation.
                new_potential_array(r - 1, z - 1) ...
                    = potential_array(floor(r / 2), floor(z / 2));
            end
        end
    
        potential_array = solve_laplace(min([d_start * double(2.0^(start - res_exp)), double(res / 64.0)]), new_potential_array, boundary_conditions, where_bounds, min(final_accuracy, 8.0));
    end
    
    imagesc(potential_array)
    set(gca,'YDir','normal')
    
    % 2 * dimensions(2) - 1 because the first cell is the zero cell
    cartesian_dimensions = ...
        [2 * dimensions(2) - 1, 2 *  dimensions(2) - 1, dimensions(1)];
    cartesian_potential_maps = ...
        zeros(cartesian_dimensions);
    
    % Grid for interpolation
    x = double(1:dimensions(2));
    y = double(1:dimensions(1));
    [X, Y] = meshgrid(x, y);
    
    electrode_map = potential_array;
    for x = (-dimensions(2) + 1):(dimensions(2) - 1)
        for y = (-dimensions(2) + 1):(dimensions(2) - 1)
            for z = double(1:dimensions(1))
                % SEE NOTES
                r = sqrt(double((abs(x) - 1)^2 + (abs(y) - 1)^2));
                % Rounding only for no interpolation case.
                r = round(r);
                if r < dimensions(2)
                    disp(z)
                    disp(r + 1.0)
                    disp(size(electrode_map))
                    cartesian_potential_maps( ...
                        dimensions(2) + x, ...
                        dimensions(2) + y, z) ...
                            = lininterp2(X, Y, electrode_map, z, r + 1.0);
                            % = potential_maps(electrode, z, r + 1);
                end
            end
        end
    end

    fly_result = test_fly(reshape(cartesian_potential_maps, [1 size(cartesian_potential_maps)]), ...
                          1.0 / 8.0, 100, 128, 128, 0, 0, 0);

    score = 0;
end

function output = test_fly(potential_array, d, x, y, z, vx, vy, vz)
    start_time =  0.0;      % us
    end_time   =  10.0;     % us
    m = 0.0006446;                % amu (e.g. 2.0 would be roughly correct for H2+)
    q = 1.0;                % atomic units
    
    maxdist =  0.00025;     % mm

    time_steps_per_us = 2;
    time_steps = round((end_time - start_time) * time_steps_per_us);
    step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
    voltages = ones([time_steps, 1]);

    tic
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    [x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
        = trajectory_integration_module(x, y, z, vx, vy, vz, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          length(electrode_names), m, q, d, maxdist, end_time);
    
    output.its = int32(its);
    output.original_length = length(ts);
    output.x_traj = x_traj(1:its);
    output.y_traj = y_traj(1:its);
    output.z_traj = z_traj(1:its);
    output.ts     =     ts(1:its);
    output.exs    =    exs(1:its);
    output.eys    =    eys(1:its);
    output.ezs    =    ezs(1:its);
    
    elapsed_time = toc;
    fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(its / elapsed_time));
end

%%
function solved = solve_laplace(d, potential_array, boundary_conditions, where_bounds, accuracy)
    dimensions = size(potential_array);
    potential_array(where_bounds) = boundary_conditions(where_bounds);
    
    change = 0.0;
    last_change = 0.0;
    i = 0;
    while (last_change - change) / (dimensions(1) * dimensions(2)) > 10.0^(-accuracy) ...
            || i < 64
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
                if ~where_bounds(r, z)
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

function [where_bounds, bounds] = generate_electrodes(vs, rs, z1s, z2s, dimensions)
    bounds = zeros(dimensions, "double");
    where_bounds = zeros(dimensions, "logical");
    for r = 1:dimensions(1)
        for z = 1:dimensions(2)
            for i = 1:length(rs)
                if r == round(rs(i) * double(dimensions(1)))
                    if z >= round(z1s(i) * double(dimensions(2))) ...
                        && z <= round(z2s(i) * double(dimensions(2)))
                        where_bounds(r, z) = true;
                        bounds(r, z) = vs(i);
                    end
                end
            end
        end
    end
end