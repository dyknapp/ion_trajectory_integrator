function [ ionStruct ] = createIonCloud( sim, radius, atomType, N, seed)
%CREATEIONCLOUD Adds a cloud of ions to the simulation. 
%  Creates N atoms of the specified atomType and adds them to the
%  simulation sim. Atoms are placed at random positions uniformly
%  distributed within a sphere of specified radius. The seed for the random
%  number generator may be specified.
% 
% SYNTAX: createIonCloud(sim, radius, atomType, N, seed)
%
% Example:
%   createIonCloud(sim, 1e-4, calciumIon, 40)

if nargin < 5
    seed = randi(100000,1);
end

rng(seed);

% Lammps does have a function that can create ions in a cloud like
% configuration; however, it requires a lattice to be declared, and is
% prone to overlapping atoms. As a result, we instead calculate positions
% in matlab and position the atoms ourselves.

x = ones(0,1);
y = x;
z = x;

for i=1:N
    %Select a random length
    d = rand * radius;
    
    %Select a random azimuthal and inclination angle
    a = pi*rand;
    b = 2*pi*rand;
    
    %Add positions to x, y, z;
    x(end+1) = d * sin(a) * cos(b);
    y(end+1) = d * sin(a) * sin(b);
    z(end+1) = d * cos(a);
end

ionStruct = placeAtoms(sim, atomType, x', y', z');
end