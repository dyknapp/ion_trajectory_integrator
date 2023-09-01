% Comments only made on code that is unique to this file.
% For general information about the code here, refer to linear_paul_trap.m
%
% dknapp, 8.30.2023: Wrote the script
%%
simion_path = "SIM001_data/001_linear_paul_trap";
electrode_names = ["quad_doubled.PA1.patxt", ...
                   "quad_doubled.PA2.patxt", ...
                   "quad_doubled.PA3.patxt", ...
                   "quad_doubled.PA4.patxt", ...
                   "quad_doubled.PA5.patxt", ...
                   "quad_doubled.PA6.patxt", ...
                   "quad_doubled.PA7.patxt", ...
                   "quad_doubled.PA8.patxt", ...
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

%% Parameter sweep values & array for results
res = 256;
amplitudes = linspace(0, 800, res);
endcap_voltages = linspace(-25, 0, res);
lifetimes = zeros(res);

%% Run parameter sweep

for a_idx = 1:res
    for e_idx = 1:res
        tic
        % Set up and run simulation
        d = 0.1;                        % mm
        start_time =  0.0;              % us
        end_time   =  10.0;             % us
        m = 2.0;                        % amu (e.g. 2.0 would be roughly correct for H2+)
        q = 1.0;                        % atomic units
        maxdist = 0.001;                % mm
        time_steps_per_us  = 1000;      % number
        RF_frequency = 15.0 * 10.0^6;   % Hz
        RF_amplitude = amplitudes(a_idx);% V
        endcap_voltage = endcap_voltages(e_idx);% V
        double_RF_electrodes = true;   % Boolean.  Should RF be applied to the center electrodes?
        
        time_steps = round((end_time - start_time) * time_steps_per_us);
        step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
        voltages = zeros([time_steps, length(electrode_names)]);
        
        % RF electrodes: 1, 2
        voltages(:, 1) = RF_amplitude * ...
            cos(2 * pi * step_times * RF_frequency / 10.0^6);
        voltages(:, 2) = RF_amplitude * ...
            cos(2 * pi * step_times * RF_frequency / 10.0^6);
        if double_RF_electrodes
            % Center electrodes: 4, 7 (pi phase shift by doing cos -> sin)
            voltages(:, 4) = RF_amplitude * ...
                sin(2 * pi * step_times * RF_frequency / 10.0^6);
            voltages(:, 7) = RF_amplitude * ...
                sin(2 * pi * step_times * RF_frequency / 10.0^6);
        end
        for electrode = [3 5 6 8]
            voltages(:, electrode) = ones([time_steps, 1], "double") * endcap_voltage;
        end
        % Endcaps: 3, 5, 6, 8
        for electrode = [3 5 6 8]
            voltages(:, electrode) = ones([time_steps, 1], "double") * endcap_voltage;
        end
        
        T = 3; % Kelvin
        maxwell = @(v) maxwell_pdf(v, m, T);
        v = general_distribution(1, 1, 10000, maxwell) * 1.0e-3; % mm / us
        theta = 2 * pi * rand();
        phi = 2 * pi * rand();
        vxx1 = v * sin(theta) * cos(phi);      % mm / us
        vyy1 = v * sin(theta) * sin(phi);      % mm / us
        vzz1 = v * cos(theta);                 % mm / us
        xx1 = d * double(dimensions(1)) / 2.0; % mm
        yy1 = d * double(dimensions(2)) / 2.0; % mm
        zz1 = d * double(dimensions(3)) / 2.0; % mm
        
        potential_maps_size = size(potential_maps);
        potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
        [x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
            = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                              potential_maps, voltages, step_times, ...
                              time_steps, dimensions, int32(is_electrode), ...
                              length(electrode_names), m, q, d, maxdist, end_time);
        elapsed = toc;
        % Save result, update status
        lifetimes(a_idx, e_idx) = ts(int32(its));
        fprintf("RF frequency:      %.1f MHz\n" + ...
                "RF amplitude:      %d V\n" + ...
                "Endcap voltage:    %d V\n" + ...
                "Ion lifetime:      %.3g us\n" + ...
                "Simulation time:   %d ms\n", ...
                RF_frequency / 1.0e+6, round(RF_amplitude), ...
                round(endcap_voltage), lifetimes(a_idx, e_idx), ...
                round(elapsed * 1.0e+3)); 
                fprintf("\n");
    end
end

%% Plot the results

% Plot the parameter sweep outcome as an image.
imagesc(amplitudes, endcap_voltages, lifetimes');
set(gca,'YDir','normal');
xlabel("RF amplitude (V)")
ylabel("Endcap voltage (V)")

%%
% To compare against Mathieu functions, we need to calculate the
%   characteristic values.  Unfortunately, this functionality is not
%   natively available in MATLAB, so I choose the somewhat awkward
%   arrangement of calling Python SciPy from MATLAB.
% You need to add the python path with an environment that has scipy
%   installed.  Use python 3.10.
% UNCOMMENT BELOW IF YOU WANT TO PLOT MATHIEU FUNCTION EIGENVALUES
% % samples = 1000;
% % pyenv('Version', 'C:\Users\h2p_l\anaconda3\envs\dknapp\python.exe', 'ExecutionMode', 'OutOfProcess');
% % pyrun('import scipy');
% % pyrun('import numpy as np');
% % qs  = double(pyrun(sprintf('qs = np.linspace(0, 2,   %d)', samples), 'qs'));
% % as  = double(pyrun(sprintf('a_s = np.linspace(0, 1, %d)', samples), 'a_s'));
% % a0s = double(pyrun('a0s = scipy.special.mathieu_a(0.0, qs)', 'a0s'));
% % b1s = double(pyrun('b1s = scipy.special.mathieu_b(1.0, qs)', 'b1s'));
% % % Plot the Mathieu functions
% % hold on
% % lw = 4;
% % plot(qs, -a0s, '--k', 'LineWidth', lw);
% % plot(qs, -b1s, '--k', 'LineWidth', lw);
% % plot(qs, b1s, '--r', 'LineWidth', lw);
% % hold off
% % legend("Mathieu -a(0, q)", "Mathieu -b(1, q)", "Mathieu  b(1, q)", "Location", "northwest")
% % xlabel('q')
% % ylabel('a')
% % title("Hyperbolic ion trap calculated stability diagram vs theoretical")

