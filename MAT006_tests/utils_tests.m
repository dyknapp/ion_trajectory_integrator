%% interpolate_monotonic_1d

n = 2^20;
m = 128;
t = sort(rand([1 n]));
f = sin(2 * pi * t) + 4 * cos(16 * pi * t) .* exp(-64 * (t - 0.5).^2);
new_t = linspace(-0.25, 1.25, m);

new_f = interpolate_monotonic_1d(f, t, new_t);

figure
hold on
plot(t, f, '-k', LineWidth=5)
plot(new_t, new_f, '.-r') 
hold off

%% interpolate_monotonic_1d speed test
% sample size upconversion
res = 16;
m_sizes = round(logspace(log10(2^5), log10(2^20), res));

times = zeros([1 res]);
for idx = 1:res
    new_t = linspace(-0.25, 1.25, m_sizes(idx));
    test_run = @() interpolate_monotonic_1d(f, t, new_t);
    times(idx) = timeit(test_run);
    fprintf("Interpolation from %d to %d samples took %.3gs\n", n, m_sizes(idx), times(idx))
end