function score = score_stack(v_repeller, ...
                             v_focus, ...
                             v_shim1, ...
                             v_lens1, ...
                             v_shim2, ...
                             v_shim3, ...
                             v_lens2, ...
                             batches, particles_per_batch, ...
                             xxs, yys, zzs, ...
                             vxxs, vyys, vzzs, ...
                             potential_maps, is_electrode, dimensions, ...
                             electrode_names)

    d = 1.0;                % mm
    start_time =  0.0;      % us
    end_time   =  1.0;     % us
    m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
    q = 1.0;                % atomic units
    maxdist =  0.1;      % mm
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
    
    % Integration
    output_points = int32(3);
    
    % Prep potential_maps
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    
    % % Initialize parallel pool
    % if isempty(gcp('nocreate'))
    %     parpool
    % end
    % % Parallel portion
    % % Distribute work for the batches
    % init_xs = reshape(xxs, [batches, particles_per_batch]);
    % init_ys = reshape(yys, [batches, particles_per_batch]);
    % init_zs = reshape(zzs, [batches, particles_per_batch]);
    % init_vxs = reshape(vxxs, [batches, particles_per_batch]);
    % init_vys = reshape(vyys, [batches, particles_per_batch]);
    % init_vzs = reshape(vzzs, [batches, particles_per_batch]);
    % % Prep array to store the results
    % result_xs = zeros([batches, particles_per_batch, output_points]);
    % result_ys = zeros([batches, particles_per_batch, output_points]);
    % result_zs = zeros([batches, particles_per_batch, output_points]);
    % % Run the computation
    % parfor i = 1:batches
    %     [txs, tys, tzs, ~, ~] = fly_ensemble(output_points, int32(particles_per_batch), ...
    %                                           init_xs(i, :),  init_ys(i, :),  init_zs(i, :), ...
    %                                           init_vxs(i, :), init_vys(i, :), init_vzs(i, :), ...
    %                                           potential_maps, voltages, step_times, ...
    %                                           int32(time_steps), dimensions, int32(is_electrode), ...
    %                                           int32(length(electrode_names)), m, q, d, ...
    %                                           maxdist, end_time);
    %     result_xs(i, :, :) = txs;
    %     result_ys(i, :, :) = tys;
    %     result_zs(i, :, :) = tzs;
    % end
    % % Reformat the results
    % xss = reshape(result_xs, [batches * particles_per_batch, output_points]);
    % yss = reshape(result_ys, [batches * particles_per_batch, output_points]);
    % zss = reshape(result_zs, [batches * particles_per_batch, output_points]);

    % Single thread run:
    [xss, yss, zss, ~, ~] = fly_ensemble( output_points, int32(n_particles), xxs, yys, zzs, vxxs, vyys, vzzs, ...
                                          potential_maps, voltages, step_times, ...
                                          int32(time_steps), dimensions, int32(is_electrode), ...
                                          int32(length(electrode_names)), m, q, d, ...
                                          maxdist, end_time);
    
    detector_radius = 20.0;
    final_rs = sqrt((xss(:, end) - mean(xxs)).^2 + (yss(:, end) - mean(yys)).^2);
    detected = (zss(:, end) <= 7.0) & (final_rs <= detector_radius);
    survival_portion = double(sum(detected)) / double(n_particles);
    score = 0.;
    % Penalty for particles not making it to the detector
    score = score + 100. * (1. - survival_portion);
    if survival_portion > (1. /double(n_particles))
        corr_matrix = corrcoef([xss(detected, end), yss(detected, end), zss(detected, end), ...
                                vxxs(detected),     vyys(detected),     vzzs(detected)]);
        for idx = 1:3
            score = score - abs(corr_matrix(idx + 3, idx));
        end
        % score = score - 3.0 * (min(mean(final_rs(detected)), detector_radius) / detector_radius);
        % % Bonus for radius of result.
        % score = score - 3. * (min(mean(final_rs(detected)), detector_radius) / detector_radius);
    end
    score = score - 0.1   * min(std(xss(detected, end)), detector_radius/2.) / (detector_radius/2.) ...
                          * min(std(xss(detected, end)), detector_radius/2.) / (detector_radius/2.);
end