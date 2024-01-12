function result = general_distribution( N, dv, vmax, pdf)
    % Generate N random doubles in the Maxwell-Boltzmann distribution.
    % Do we need to take into account rovibrational dof's?
    % INPUT:
    %       N       : How many numbers to output
    %       dv      : Our result is actually discretized.  How large should the
    %                   bins be?
    %       vmax    : Where to cut off the distribution

    dv = double(dv);
    vmax = double(vmax);
    n_bins = round(vmax / dv);
    v_bins = double(1:n_bins) * (vmax / (n_bins + 1));
    
    % Iterate over v_bins and "integrate" to get CDF
    cdf = zeros([length(v_bins), 1]);
    cdf(1) = pdf(v_bins(1));
    for v_idx = 2:length(v_bins)
        cdf(v_idx) = cdf(v_idx - 1) + dv * pdf(v_bins(v_idx));
    end
    cdf = (cdf - cdf(1));
    cdf = cdf  / cdf(end);

    % Take uniform random variables and interpolate where they'd end up in
    % the cdf.
    uniform_randoms = rand(N, 1) * (1 - 1.0e-15);
    result = zeros([N, 1]);

    % Since interp1 doesn't tolerate repeats in cdf
    [cdf, idcs] = unique(cdf,'first');
    v_bins = v_bins(idcs);

    for idx = 1:N
        ar_idx = find((cdf - uniform_randoms(idx)) > 0, 1, "first");
        br_idx = ar_idx - 1;
        if br_idx == 0
            result(idx) = uniform_randoms(idx) * (v_bins(2) - v_bins(1)) / cdf(1);
        else
            result(idx) = v_bins(br_idx) +  (v_bins(ar_idx) - v_bins(br_idx))...
                                            * ((uniform_randoms(idx) - cdf(ar_idx)) ...
                                            / (cdf(ar_idx) - cdf(br_idx)));
        end
    end
end
