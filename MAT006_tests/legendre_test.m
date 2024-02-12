res = 16384;
xs = linspace(-1, 1, res);
ys = zeros(size(xs));

l = 64;
for idx = 1:res
    ys(idx) = legendre_f(int32(l), int32(0), double(xs(idx)));
end
plot(xs, ys, '.');

hold on
matlab_ys = legendre(l, xs);
plot(xs, matlab_ys(1,:), '-')
hold off

xlim([-1.1, 1.1])
ylim([-1.1, 1.1])

fprintf('RMS error %.3g\n', mean((ys - matlab_ys(1,:)).^2));