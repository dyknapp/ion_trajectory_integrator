res = 1024;
n_terms = 128;

% Boundary conditions:
a_cond =  0.0;
b_cond =  3.0;

xs = linspace(0.0, 2 * pi, res);

% % Plot the lifting function
% plot(xs, lift(potential_is_one_between, xs))

%% Calculate the basis functions ahead of time

grad_xs = gradient(xs);
basis_calc = zeros([n_terms, res]);
grad_basis_calc = zeros([n_terms, res]);
parfor idx = 1:n_terms
    basis_calc(idx, :) = basis(idx, xs);
    grad_basis_calc(idx, :) = basis_calc(idx, :) ./ grad_xs;
end

%% Construct the A matrix

A = zeros(n_terms);
% Symmetric, so fill upper triangle and then symmetrize
for idx1 = 1:n_terms
    for idx2 = 1:n_terms
        A(idx1, idx2) = ...
            sum(grad_xs.^2 .* grad_basis_calc(idx1, :) ...
            .* grad_basis_calc(idx2, :), "all");
    end
end

%% Construct the b vector

c = zeros([n_terms 1]);
lift_derivative = gradient(lift(a_cond, b_cond, xs)) ./ grad_xs;
for idx = 1:n_terms
    c(idx) = ...
        sum(grad_xs.^2 .* grad_basis_calc(idx, :) ...
            .* lift(a_cond, b_cond, xs), "all");
        % .* lift_derivative, "all");
end

b = zeros([n_terms 1]);
spike_calc = spike(xs);
for idx = 1:n_terms
    b(idx) = ...
        sum(grad_xs.^2 .* basis_calc(idx, :) ...
        .* spike_calc, "all");
end

b = b - A * c;

%% Solve for the coefficients

coeffs = linsolve(A, b);

%% Reconstruct the potential

potential = zeros([1 res]);
for idx = 1:n_terms
    potential = potential + coeffs(idx) * basis_calc(idx, :);
end

plot(xs, potential);

%%

function ys = lift(a_cond, b_cond, xs)
    ys = a_cond * xs + (b_cond * (1.0 - max(xs)));
end

function ys = basis(n, xs)
    ys = besselj(n, xs);
end

function rho = spike(xs)
    rho = exp(-100 * (xs - 0.25).^2);
end