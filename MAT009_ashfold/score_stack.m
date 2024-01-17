function score = score_stack(v_repeller, ...
                             v_focus, ...
                             v_shim1, ...
                             v_lens1, ...
                             v_shim2, ...
                             v_shim3, ...
                             v_lens2, ...
                             batches, particles_per_batch, ...
                             vtransverse, ...
                             xxs, yys, zzs, ...
                             potential_maps, is_electrode, dimensions, ...
                             electrode_names)

    d = 1.0;                % mm
    start_time =  0.0;      % us
    end_time   =  1.0;     % us
    m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
    q = 1.0;                % atomic units
    maxdist =  0.01;      % mm
    n_particles = batches * particles_per_batch;  % 2^10;

    time_steps_per_us = 1000;
    time_steps = round((end_time - start_time) * time_steps_per_us);
    step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
    voltages = zeros([time_steps, length(electrode_names)]);
    turn_on_time = 0.1;

    % Always have the detector at -2.5kV
    voltages = set_voltage_at_time(8, -2500.0,    0.0, step_times, voltages);

    voltages = set_voltage_at_time(1, v_repeller, turn_on_time, step_times, voltages);
	voltages = set_voltage_at_time(2, v_focus,    0.0, step_times, voltages);
	voltages = set_voltage_at_time(3, v_shim1,    0.0, step_times, voltages);
	voltages = set_voltage_at_time(4, v_lens1,    0.0, step_times, voltages);
	voltages = set_voltage_at_time(5, v_shim2,    0.0, step_times, voltages);
	voltages = set_voltage_at_time(6, v_shim3,    0.0, step_times, voltages);
	voltages = set_voltage_at_time(7, v_lens2,    0.0, step_times, voltages);

    % Initialize
    vxx1 = vtransverse / sqrt(2.);
    vyy1 = vtransverse / sqrt(2.);
    vzz1 = 0.0;
    vxxs = vxx1 * ones([n_particles 1]);
    vyys = vyy1 * ones([n_particles 1]);
    vzzs = vzz1 * ones([n_particles 1]);
    
    % Integration
    output_points = int32(3);
    
    % Initialize parallel pool
    if isempty(gcp('nocreate'))
        parpool
    end
    % Prep potential_maps
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    
    % Parallel portion
    % Distribute work for the batches
    init_xs = reshape(xxs, [batches, particles_per_batch]);
    init_ys = reshape(yys, [batches, particles_per_batch]);
    init_zs = reshape(zzs, [batches, particles_per_batch]);
    init_vxs = reshape(vxxs, [batches, particles_per_batch]);
    init_vys = reshape(vyys, [batches, particles_per_batch]);
    init_vzs = reshape(vzzs, [batches, particles_per_batch]);
    % Prep array to store the results
    result_xs = zeros([batches, particles_per_batch, output_points]);
    result_ys = zeros([batches, particles_per_batch, output_points]);
    result_zs = zeros([batches, particles_per_batch, output_points]);
    % Run the computation
    parfor i = 1:batches
        [txs, tys, tzs, ~, ~] = fly_ensemble(output_points, int32(particles_per_batch), ...
                                              init_xs(i, :),  init_ys(i, :),  init_zs(i, :), ...
                                              init_vxs(i, :), init_vys(i, :), init_vzs(i, :), ...
                                              potential_maps, voltages, step_times, ...
                                              int32(time_steps), dimensions, int32(is_electrode), ...
                                              int32(length(electrode_names)), m, q, d, ...
                                              maxdist, end_time);
        result_xs(i, :, :) = txs;
        result_ys(i, :, :) = tys;
        result_zs(i, :, :) = tzs;
    end
    % Reformat the results
    xss = reshape(result_xs, [batches * particles_per_batch, output_points]);
    yss = reshape(result_ys, [batches * particles_per_batch, output_points]);
    zss = reshape(result_zs, [batches * particles_per_batch, output_points]);
    % [xss, yss, zss, ~, ~] = fly_ensemble( output_points, int32(n_particles), xxs, yys, zzs, vxxs, vyys, vzzs, ...
    %                                       potential_maps, voltages, step_times, ...
    %                                       int32(time_steps), dimensions, int32(is_electrode), ...
    %                                       int32(length(electrode_names)), m, q, d, ...
    %                                       maxdist, end_time);
    
    detected = zss(:, end) <= 7.0;
    survival_portion = double(sum(detected)) / double(n_particles);
    score = 0;
    if survival_portion > 0
        mean_r = mean(sqrt((xss(detected, end) - mean(xxs(detected, end))).^2 ...
            + (yss(detected, end) - mean(yys(detected, end))).^2));
        % score = score + sqrt(std(xss(detected, end))^2. + std(yss(detected, end))^2.) / mean_r;
        score = score - mean_r;
    end
    % Penalty for particles not making it to the detector
    score = score + (1 - survival_portion) * (10000. / 3.);
end