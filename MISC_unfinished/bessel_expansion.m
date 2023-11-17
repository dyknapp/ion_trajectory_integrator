res = 1024;
r_max = 3.0;
n_max = 10;
p = 1;

f = @(r) r.^2 .* (1 - r / r_max);

rs = linspace(0, r_max, res);

figure;
plot(rs, f(rs), '-k', 'LineWidth', 5);

%% Expansion coeffs

% initialize zeros of Bessel J function
bessel_zeros = besselzero(p, n_max, 1);

coeffs = zeros([1 n_max]);
for n = 1:n_max
    integrand = @(r) (r .* f(r) .* besselj(p, bessel_zeros(n) .* r ./ r_max));
    coeffs(n) = integral(integrand, 0, r_max);
    coeffs(n) = coeffs(n) * 2 / (r_max * besselj(p + 1, bessel_zeros(n)))^2;
end

%% Plot expansion

expansion = zeros([1 res]);
for n = 1:n_max
    expansion = expansion + coeffs(n) * besselj(p, bessel_zeros(n) .* rs ./ r_max);
end

hold on
plot(rs, expansion, '-g', 'LineWidth', 2);
hold off

title(sprintf('mean squared error = %.3g', sum((expansion - f(rs)).^2) / r_max))