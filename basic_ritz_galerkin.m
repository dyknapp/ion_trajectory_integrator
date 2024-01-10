res = 128;
extent = 16.0;

xs = linspace(-extent, extent, res);
ys = linspace(-extent, extent, res);

% % The charge distribution:
% plot_img(xs, ys, f(xs, ys));

% % The basis functions
% m = -0.5; n = -0.5;
% plot_img(xs, ys, basis(m, n, extent, xs, ys))

% n_terms > 1, must be an integer
index_terms = 32;

m = (0:index_terms - 1);
n = m;

ms = zeros([index_terms^2 1]);
ns = ms;
count = 0;
for idx1 = 1:index_terms
    for idx2 = 1:index_terms
        count = count + 1;
        ms(count) = m(idx1);
        ns(count) = n(idx2);
    end
end

n_terms = index_terms^2;

%%  Construct matrix A

dxdy = (2 * extent / (res - 1))^2;

tic
A = zeros(n_terms);
parfor idx1 = 1:n_terms
    for idx2 = 1:n_terms
        m1 = ms(idx1);
        n1 = ns(idx1);
        m2 = ms(idx2);
        n2 = ns(idx2);
        A(idx1, idx2) = dxdy ...
            * sum(grad_basis(m1, n1, extent, xs, ys) ...
               .* grad_basis(m2, n2, extent, xs, ys), "all");
    end
end
toc

%% Construct vector b

tic
b = zeros([n_terms 1]);
parfor idx1 = 1:n_terms
    m1 = ms(idx1);
    n1 = ns(idx1);
    b(idx1) = dxdy ...
        * sum(basis(m1, n1, extent, xs, ys) ...
               .* f(xs, ys), "all");
end
toc

% %% extra boundary conditions
% 
% corner_zero = 4.0;
% c = zeros([n_terms 1]);
% parfor idx1 = 1:n_terms
%     m1 = ms(idx1);
%     n1 = ns(idx1);
%     b(idx1) = dxdy ...
%         * sum(basis(m1, n1, extent, xs, ys) ...
%               .* double(xs < corner_zero) .* double(ys < corner_zero), "all");
% end
% b = b - A*c;

%% Solve for coefficients

tic
coeffs = linsolve(A, b);
toc

%% Calculate the potential

tic
potential = zeros(res);
parfor idx1 = 1:n_terms
    m1 = ms(idx1);
    n1 = ns(idx1);
    potential = potential ...
        + coeffs(idx1) * basis(m1, n1, extent, xs, ys);
end
toc
plot_img(xs, ys, potential);

%%

% imagesc(A); axis image;

%% Relevant functions
% y is transposed so that we get a matrix out for 1D x and y
% e.g.
% 
% >> [10 20 30] .* [1 2 3]'
% 
% ans =
% 
%     10    20    30
%     20    40    60
%     30    60    90


function density = f(x, y)
    density = 1.0 * double((x.^2 + (y.^2)') <= 8);
end

function potential = basis(m, n, extent, x, y)
    % m,n = 0, 1, 2, ...
    km = pi * double(2 * m + 1) / (2.0 * extent);
    kn = pi * double(2 * n + 1) / (2.0 * extent);
    potential = cos(km * x) .* cos(kn * y)';
end

function grad_potential = grad_basis(m, n, extent, x, y)
    % m,n = 0, 1, 2, ...
    km = pi * double(2 * m + 1) / (2.0 * extent);
    kn = pi * double(2 * n + 1) / (2.0 * extent);
    grad_potential = [(km * -sin(km * x)) .* cos(kn * y)', ...
                      cos(km * x) .* (kn * -sin(kn * y)')];
end