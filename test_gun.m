vs = [-15, 100, -5, 0];
r1s = [0.50, 0.30, 0.50, 0.50];
r2s = [1.0 1.0 1.0 1.0];
z1s = [0.01, 0.1, 0.2, 0.3];
z2s = [0.09, 0.19, 0.39, 0.4];

tic
test(100, 300, vs, r1s, r2s, z1s, z2s)
toc

%%

function score = test(res_r, res_z, vs, r1s, r2s, z1s, z2s)
    dimensions = [res_r+2, res_z+2];
    potential_array = zeros(dimensions);
    boundary_conditions = zeros(dimensions);
    guess = zeros(dimensions);
    
    [bc_mask, bcs] = generate_electrodes(vs, r1s, r2s, z1s, z2s, dimensions);
    guess(bc_mask == 1.0) = bcs(bc_mask == 1.0);


    % Solve Laplace
    refined = refined_laplace(guess, 1.0 - bc_mask, 1.0e-6, 1024, dimensions(1), dimensions(2), 4);
    % Cut away the "ghost cells" along the border of the mesh
    refined = refined(2:end-1, 2:end-1);

    imagesc(refined)
    set(gca,'YDir','normal')
    
    % 2 * dimensions(2) - 1 because the first cell is the zero cell
    cartesian_dimensions = ...
        [2 * dimensions(2) - 1, 2 *  dimensions(2) - 1, dimensions(1)];
    cartesian_potential_maps = ...
        zeros(cartesian_dimensions);

    % Grid for interpolation
    X = double(1:(dimensions(1) - 2));
    Y = double(1:(dimensions(2) - 2));
    
    cartesian = zeros([2*res_r + 1, 2*res_r + 1, res_z]);
    for x_grid = 1:(2*res_r + 1)
        for y_grid = 1:(2*res_r + 1)
            for z_grid = 1:res_z
                x_coord = double(x_grid) - double(res_r + 1);
                y_coord = double(y_grid) - double(res_r + 1);
                r_coord = sqrt(x_coord^2 + y_coord^2);
                r_grid = r_coord + 1.0;
                if r_grid < res_r
                    value = lininterp2(X, Y, refined, r_grid, double(z_grid));
                    cartesian(x_grid, y_grid, z_grid) = value;
                end
            end
        end
    end

    fly_result = test_fly(reshape(cartesian, [1 size(cartesian)]), zeros([2*res_r + 1, 2*res_r + 1, res_z]),...
                          0.1, 100, double(res_r), double(res_r), 0, 0, 0);

    score = 0;
end

function output = test_fly(potential_maps, is_electrode, d, x, y, z, vx, vy, vz)
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
    dimensions = potential_maps_size(2:4);
    [x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
        = trajectory_integration_module(x, y, z, vx, vy, vz, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          1, m, q, d, maxdist, end_time);
    
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
function [where_bounds, bounds] = generate_electrodes(vs, r1s, r2s, z1s, z2s, dimensions)
    bounds = zeros(dimensions, "double");
    where_bounds = zeros(dimensions, "logical");
    for r = 1:dimensions(1)
        for z = 1:dimensions(2)
            for i = 1:length(vs)
                if r >= round(r1s(i) * double(dimensions(1)))...
                   && r <= round(r2s(i) * double(dimensions(2)))
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