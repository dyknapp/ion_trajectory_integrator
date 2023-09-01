function result = maxwell_distribution(N, m, T, dv, vmax)
    % Generate N random doubles in the Maxwell-Boltzmann distribution.
    % Do we need to take into account rovibrational dof's?
    % INPUT:
    %       N       : How many numbers to output
    %       m       : Particle mass in amu
    %       T       : Temperature in Kelvin
    %       dv      : Our result is actually discretized.  How large should the
    %                   bins be?
    %       vmax    : Where to cut off the distribution
    
    m = double(m);
    T = double(T);
    dv = double(dv);
    vmax = double(vmax);
    n_bins = round(vmax / dv);
    v_bins = double(1:n_bins) * (vmax / (n_bins + 1));
    
    % Iterate over v_bins and "integrate" to get CDF
    cdf = zeros([length(v_bins), 1]);
    cdf(1) = maxwell_pdf(v_bins(1), m, T);
    for v_idx = 2:length(v_bins)
        cdf(v_idx) = cdf(v_idx - 1) + dv * maxwell_pdf(v_bins(v_idx), m, T);
    end
    
    % Take uniform random variables and interpolate where they'd end up in
    % the cdf.
    uniform_randoms = rand(N, 1) / cdf(end);
    result = zeros([N, 1]);

    % Since interp1 doesn't tolerate repeats in cdf
    [cdf, idcs] = unique(cdf,'first');
    v_bins = v_bins(idcs);

    for idx = 1:N
        result(idx) = interp1(cdf, v_bins, uniform_randoms(idx));
    end
end