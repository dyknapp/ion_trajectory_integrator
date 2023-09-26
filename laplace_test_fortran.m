res_r = 100;
res_z = 200;
dimensions = [res_r+2, res_z+2];
potential_array = zeros(dimensions);
boundary_conditions = zeros(dimensions);
guess = zeros(dimensions);

% JUST THIS ONE IS ONES:
bc_mask = ones(dimensions);

% Boundary conditions
for r = 2:dimensions(1) - 1
    for z = 2:dimensions(2) - 1
        if r >= 41 && r <= 61
            boundary_conditions(r, z) = 1.0;
            bc_mask(r, z) = 0.0;
        end
    end
end

% Apply initial conditions
guess(bc_mask == 0.0) = boundary_conditions(bc_mask == 0.0);

tic
[its, refined] = iterate_laplace(guess, bc_mask, 1.0e-4, 1024, dimensions(1), dimensions(2));
disp(its);
toc
imagesc(refined(2:end-1, 2:end-1));
set(gca,'YDir','normal')
colorbar()

%%

tic
refined = refined_laplace(guess, bc_mask, 1.0e-7, 2^14, dimensions(1), dimensions(2), 3);
% Cut away the "ghost cells" along the border of the mesh
refined = refined(2:end-1, 2:end-1);
toc
imagesc(refined);
set(gca,'YDir','normal')
colorbar()






