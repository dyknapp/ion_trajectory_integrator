function [pot_array, is_electrode, dimensions] = ...
    readFile(potpathlist, startline)
% function to read the potentials from a patxt file to matlab
%
% mbroerse, 16.6.2023
% dknapp,   14.8.2023
%
% INPUT:
%       potpathlist ... vector of paths to .patxt files
%       startline   ... integer that represents the first line of data

% OUTPUT:
%       pot_array       ... array of 3D potentials in V for .patxt files
%                           generated from .pa0 files and in 10^-4 V for
%                           potentials generated from .paX files. this array
%                           has the same length as potpathlist.
%       is_electrode    ... 3D array for wich the presence of an electrode
%                           is indicated with 1, otherwise 0.
%       dimensions      ... vector of length 3 which indicate the size of
%                           the abovementioned arrays in the x, y, z
%                           direction.
%                           
        

    % get first array of potentials in the format: 'x y z is_electrode
    % potential', start with the header
    fid = fopen(potpathlist(1));
    scan = textscan(fid, '%*s %d', 3, 'Headerlines', 6);
    dimensions = scan{1}.'; % can be shorter?

    % dimensions found, initialize pot_array
    pot_array = zeros([length(potpathlist), prod(dimensions)], "double");
    
    % scan the first pot_file to get is_electrode the first pot_array 
    frewind(fid);
    scan = textscan(fid, '%*d %*d %*d %d %f', prod(dimensions), 'HeaderLines', startline-1);
    [is_electrode, pot_array(1,:)] = scan{1,:};
    fclose(fid);
    
    % get the rest of the pot_arrays 
     for i = 2:length(potpathlist) 
        fid = fopen(potpathlist(i));
        scan = textscan(fid, ' %*d %*d %*d %*d %f', prod(dimensions), 'HeaderLines', startline-1); % could also throw away string
        pot_array(i,:) = scan{1,:};
        fclose(fid);
     end
    
    % reshape is_electrode for coord system, pot_array will be reshaped in
    % the getPotential function
    is_electrode = reshape(is_electrode, dimensions);
end

