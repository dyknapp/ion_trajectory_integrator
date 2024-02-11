name = 'paul_trap' ;
% name = 'Stanford_Bunny';

%% OPENSCAD rendering + stl->msh conversion
system("openscad.com -o paul_trap.stl paul_trap_parametric.scad -D electrode_number=1 -D param1=10 -D param2=1 -D param3=1");
geom = loadstl([name '.stl']);
savemsh(['files/' name '.msh'], geom);
clear geom

%% JIGSAW library config

%------------------------------------ setup files for JIGSAW
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

%------------------------------------ read GEOM. for display
geom = loadmsh (opts.geom_file);
figure ; drawmesh(geom);
view(-30,30); axis image;
title('INPUT GEOMETRY');

%% Refine geometry

opts.verbosity = 1;
opts.mesh_kern = 'delfront';
% opts.mesh_kern = 'delaunay';
opts.mesh_dims = 2 ;
opts.hfun_hmax = 0.01 ;
opts.geom_feat = true ;             % do sharp feat.'s
opts.mesh_top1 = true ;

mesh = jigsaw  (opts) ;

figure ; drawmesh(mesh);
view(-30,30); axis image;
title('JIGSAW (KERN=delfront)') ;

%% Triangles

triangles = mesh.tria3.index(:, 1:3);
points    = mesh.point.coord(:, 1:3);

triangles_size = size(triangles);
areas = zeros([1 triangles_size(1)]);
% tic
for idx = 1:triangles_size(1)
    idcs = triangles(idx, :);
    tri_pts = points(idcs, :);
    areas(idx) = ...
        0.5 * norm(cross(tri_pts(2,:)-tri_pts(1,:), tri_pts(3,:)-tri_pts(1,:)));
end
% toc

