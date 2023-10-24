function lifetime = try_real_trap_aq(RF_frequency, RF_amplitude, offset_volts, endcap_volts, dimensions, is_electrode, potential_maps)
    electrode_names = ["einzel3" "einzel2" "einzel1" "endcap1" "endcap2" "rfelecs"];
    
    % RF_frequency = 13.2e+6;
    % RF_amplitude = 270.0;
    % offset_volts = 0.0;
    % endcap_volts = 5.0;
    
    m = 1.00727647 + 2.01410178 + 0.00054858;
    q = 1.0;
    
    end_time = 100.0;
    maxdist  = 2.5e-3;
    
    d = 0.125;
    
    xx1 =      d * double(dimensions(1) + 1) / 2.0;
    yy1 =      d * double(dimensions(2) + 1) / 2.0;
    zz1 =      d * 400.0;
    
    vxx1 = normrnd(0.0, 0.05);
    vyy1 = normrnd(0.0, 0.05);
    vzz1 = normrnd(0.0, 10.0);
    
    % RF and endcap electrode voltages.
    time_steps_per_us = 100;
    time_steps = round((end_time - 0.0) * time_steps_per_us);
    step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
    voltages = zeros([time_steps, length(electrode_names)]);
    
    rf_electrode = 6;
    endcaps_in = 4;
    endcaps_out = 5;
    
    voltages(:, rf_electrode) = RF_amplitude * ...
        cos(2 * pi * step_times * RF_frequency / 10.0^6);
    voltages(:, rf_electrode) = voltages(:, rf_electrode) + offset_volts;
    voltages(:, endcaps_in) = endcap_volts * ones([time_steps 1]);
    voltages(:, endcaps_out) = endcap_volts * ones([time_steps 1]);
    
    % Simulation
    fprintf("RF frequency = %.1f MHz\n" + ...
            "RF amplitude = %.1f V\n" + ...
            "DC offset    = %.1f V\n", ...
            RF_frequency / 1.0e+6, ...
            RF_amplitude, ...
            offset_volts);
    tic
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    [~, ~, ~, ts, ~, ~, ~, its] ...
        = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          length(electrode_names), m, q, d, maxdist, end_time);
    elapsed_time = toc;
    
    its = int32(its);
    
    fprintf("Simulation took %.3gs (%d it/s)\n", ...
        elapsed_time, round(its / elapsed_time));
    lifetime = ts(its);
end 




