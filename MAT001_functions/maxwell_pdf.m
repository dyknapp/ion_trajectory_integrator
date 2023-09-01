function p = maxwell_pdf(v, m, T)
    % Constants and unit conversions
    m = (1.660539067e-27) * m;  % amu to kg
    k = 1.380649e-23;           % J / K

    p = (m / (2 * pi * k * T))^(1.5) * (4 * pi * v.^2) ...
        .* exp(- m * v.^2 / (2 * k * T));
end