function output = ...
    try_mathieu_aq(a, q, m, RF_frequency, T, end_time, maxdist, ...
    initialize, xx1, yy1, zz1, vxx1, vyy1, vzz1)
    
    % a and q parameters
    r0 = 1.0e-3;
    e  = 1.602176634e-19;
    
    omega = RF_frequency * (2 * pi);
    endcap_voltage =  a * (((1.660539067e-27) * m) * omega^2 * r0^2) / (2 * e);
    RF_amplitude   = -q * (((1.660539067e-27) * m) * omega^2 * r0^2) / (    e);

    fprintf("RF frequency = %.1f MHz\n" + ...
            "RF amplitude = %.1f V\n" + ...
            "Endcap volts = %.1f V\n" + ...
            "a            = %.3g\n" + ...
            "q            = %.3g\n", ...
            RF_frequency / 1.0e+6, ...
            RF_amplitude, ...
            endcap_voltage, ...
            a, q);

    % Variables to set beforehand
    simion_path = "SIM001_data/001_linear_paul_trap";
    electrode_names = ["hyperbolic_trap.pa1.patxt", ...
                       "hyperbolic_trap.pa2.patxt"
                      ];
    start_line = 22;
    loadanyways = false;
    if ~exist('dimensions', 'var') || loadanyways
        addpath(simion_path)
        [potential_maps, is_electrode, dimensions] = ...
            readFile(electrode_names, start_line);
        potential_maps = potential_maps / 10000.0;
        is_electrode = logical(is_electrode);
    end
    
    d = 0.1;                % mm
    start_time =  0.0;      % us
    m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
    q = 1.0;                % atomic units
    
    
    % RF and endcap electrode voltages.
    time_steps_per_us = 1000;
    time_steps = round((end_time - start_time) * time_steps_per_us);
    step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
    voltages = zeros([time_steps, length(electrode_names)]);
    
    
    % RF electrode: 2
    rf_electrode = 2;
    voltages(:, rf_electrode) = RF_amplitude * ...
        cos(2 * pi * step_times * RF_frequency / 10.0^6);
    voltages(:, rf_electrode) = voltages(:, rf_electrode) + endcap_voltage;
    
    if ~initialize
        xx1 = d * double(dimensions(1)) / 2.0; % mm
        yy1 = d * double(dimensions(2)) / 2.0; % mm
        zz1 = d * double(dimensions(3)) / 2.0; % mm
    
        maxwell = @(v) maxwell_pdf(v, m, T);
        v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3;
        theta = 2 * pi * rand();
        phi = 2 * pi * rand();
        vxx1 = v * sin(theta) * cos(phi);      % mm / us
        vyy1 = v * sin(theta) * sin(phi);      % mm / us
        vzz1 = v * cos(theta);                 % mm / us
    end

    % Simulation
    tic
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    [x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
        = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          length(electrode_names), m, q, d, maxdist, end_time);
    output.elapsed_time = toc;

    output.its = int32(its);
    original_length = length(ts);
    output.x_traj = x_traj(1:its);
    output.y_traj = y_traj(1:its);
    output.z_traj = z_traj(1:its);
    output.ts     =     ts(1:its);
    output.exs    =    exs(1:its);
    output.eys    =    eys(1:its);
    output.ezs    =    ezs(1:its);

    fprintf("Simulation took %.3gs (%d it/s)\n", ...
        output.elapsed_time, round(its / output.elapsed_time));
    output.lifetime = output.ts(end);
end