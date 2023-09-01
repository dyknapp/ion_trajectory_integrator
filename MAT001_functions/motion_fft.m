function motion_fft(ts, motion_data, RF_frequency, RF_amplitude, motion_name)

% Remove constant offset
    clean_traj = motion_data - mean(motion_data);
    % Resample so that points are equally spaced
    new_samples = linspace(ts(1), ts(end), ...
                           round(length(clean_traj) / 2) * 2);
    dt = new_samples(2) - new_samples(1);
    clean_traj = interp1(ts, clean_traj, new_samples);
    
    L = double(length(clean_traj));
    P2 = abs(fft(clean_traj) / L);
    P1 = P2(1:L / 2 + 1);
    P1(2:end - 1) = 2 * P1(2:end - 1);
    
    f = (1.0 / dt) * (0:(L/2))/L;
    figure
    plot(f,P1, '-')
    xlim([0, RF_frequency * 2.0e-6])
    xlabel('Frequency (MHz)')
    ylabel('Amplitude (mm)')
    title(sprintf('Frequency components of ion motion %s component', motion_name))
    xline(RF_frequency * 1.0e-6, '--')
    max_index = find(P1 == max(P1), 1);
    max_freq = f(max_index) * 1.0e+6;
    xline(RF_frequency * 1.0e-6 - max_freq * 1.0e-6, '--')
    xline(RF_frequency * 1.0e-6 + max_freq * 1.0e-6, '--')
    legend('FFT', 'RF Frequency', 'RF - Secular', 'RF + Secular')
    
    % Find peak frequency component -> calculate "harmonic" potential shape
    fprintf('Secular x frequency:                 %.3g MHz\n', max_freq / 1.0e+6);
    
    % r_0 calculation -> see notes (???)
    r0 = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27^2 * RF_frequency^2 * max_freq^2))^0.25;
    fprintf("r_0 predicted from FFT result is:    " + ...
            "%.3g mm\n", 1.0e+3 * r0);
    
    % Trap depth
    trap_depth = ((1.602e-19^2 * RF_amplitude^2) / (4.0 * 1.6605e-27 * RF_frequency^2 * r0^2)) / 1.602e-19;
    fprintf("Trap depth is:                       " + ...
            "%.3g eV (%.3g K)\n", trap_depth, trap_depth / 8.617e-5);

end