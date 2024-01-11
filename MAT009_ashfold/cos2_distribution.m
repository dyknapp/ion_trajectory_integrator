cos2 = @(theta) (cos(theta).^2.) / (pi);

integral(cos2, 0, 2*pi)

figure
histogram(general_distribution( 100000, 0.001, 2*pi, cos2))
