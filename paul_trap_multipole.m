n_electrodes = 4;
electrode_voltages = [1 -1 1 -1];
electrodes_file = 'paul_trap_parametric';
scad_file       = 'paul_trap_parametric.scad';

%% Generate electrodes from OpenSCAD

opts.verbosity = 1;
opts.mesh_kern = 'delfront';
% opts.mesh_kern = 'delaunay';
opts.mesh_dims = 2 ;
opts.hfun_hmax = 0.05 ;
opts.geom_feat = true ;             % do sharp feat.'s
opts.mesh_top1 = true ;

meshes = cell([1 n_electrodes]);
for idx = 1:n_electrodes
    name = sprintf('%s_%d', electrodes_file, idx);
    system(sprintf('openscad.com -o %s.stl %s -D electrode_number=%d -D param1=10 -D param2=1 -D param3=1', name, scad_file, idx));
    geom = loadstl([name '.stl']);
    savemsh(['files/' name '.msh'], geom);
    rootpath = fileparts( ...
        mfilename( 'fullpath' )) ;
    opts.geom_file = ...                % domain file
        fullfile(rootpath,...
        'files',[name,'.msh']) ;
    opts.jcfg_file = ...                % config file
        fullfile(rootpath,...
        'cache',[name,'.jig']) ;
    opts.mesh_file = ...                % output file
        fullfile(rootpath,...
        'cache',[name,'.msh']) ;
    initjig ;                           % init jigsaw

    meshes{idx} = jigsaw  (opts) ;
end

%% Consolidate the meshes

triangles = meshes{1}.tria3.index(:, 1:3);
points    = meshes{1}.point.coord(:, 1:3);
for idx = 2:n_electrodes
    triangles_size = size(triangles);
    electrode_border(idx - 1) = triangles_size(1);
    points_size = size(points);
    triangles = [triangles; (meshes{idx}.tria3.index(:, 1:3) + points_size(1))];
    points    = [points;     meshes{idx}.point.coord(:, 1:3)];
end
triangles_size = size(triangles);
electrode_border(idx) = triangles_size(1);

%% Set up the linear problem

n_triangles = size(triangles);
n_triangles = n_triangles(1);

voltages = [1 0 0 0];
voltage_vector = zeros([n_triangles 1]);
so_far = 0;
for idx = 1:n_electrodes
    voltage_vector((so_far+1):electrode_border(idx)) = voltages(idx);
    so_far = electrode_border(idx);
end

% Centroids of triangles
coords = zeros([n_triangles, 3]);
for idx = 1:n_triangles
    coords(idx, :) = sum(points(triangles(idx, :), :)) / 3.;
end

% Triangle areas
triangles_size = size(triangles);
areas = zeros([1 triangles_size(1)]);
for idx = 1:triangles_size(1)
    idcs = triangles(idx, :);
    tri_pts = points(idcs, :);
    areas(idx) = ...
        0.5 * norm(cross(tri_pts(2,:)-tri_pts(1,:), tri_pts(3,:)-tri_pts(1,:)));
end

% Matrix initialization
tic;
interaction_matrix = zeros(n_triangles);
for idx = 1:n_triangles
    line_length = fprintf("Matrix initialization progress: %.1f%%", 100.0 * double(idx) / double(n_triangles));
    for jdx = (idx+1):n_triangles
        dist = norm(coords(idx, :) - coords(jdx, :));
        value = areas(idx) * areas(jdx) / dist;
        interaction_matrix(idx, jdx) = value;
    end
    interaction_matrix(idx, idx) = 0.25 * areas(idx);
    fprintf(repmat('\b',1,line_length));
end
elapsed = toc;
fprintf("\nTime taken for matrix initialization: %.3g\n", elapsed);
interaction_matrix = interaction_matrix + interaction_matrix';
fprintf("Interaction matrix symmetric: %d\n", ...
    issymmetric(interaction_matrix));

%% Solve

tic
solution = mldivide(interaction_matrix, voltage_vector);
elapsed = toc;
fprintf("\nTime taken for matrix solution: %.3g\n", elapsed);

% interaction_matrix_gpu = gpuArray(interaction_matrix);
% solution = mldivide(interaction_matrix_gpu, voltage_vector);

%% Plot

to_plot = 1:electrode_border(1);
scatter3(coords(to_plot, 1), coords(to_plot, 2), coords(to_plot, 3), ones(size(solution(to_plot))), solution(to_plot), 'O')
axis image