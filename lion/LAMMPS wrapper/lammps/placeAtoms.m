function [ atoms ] = placeAtoms(sim, atomType, x, y, z )
%PLACEATOMS Places the given atomType at the x, y, z coordinates specified.
% Vector quantities x, y, z must be specified as column vectors.
% 
% Syntax: placeAtoms( sim, atomType, x, y, z )
%
% See Also: AddAtomType

ip = inputParser;
ip.addRequired('sim', @(sim) isa(sim, 'LAMMPSSimulation') && numel(sim) == 1);
ip.addRequired('atomType', @(type) isa(type, 'AtomType') && numel(atomType) == 1);
ip.addRequired('x', @(arg) iscolumn(arg));
ip.addRequired('y', @(arg) iscolumn(arg));
ip.addRequired('z', @(arg) iscolumn(arg));
ip.parse(sim, atomType, x, y, z);

if numel(x) ~= numel(y) || numel(x) ~= numel(z)
   error('x, y, z must be column vectors of the same length');
end

atoms = arrayfun(@(x,y,z) AtomPlacement(sim, x, y, z, atomType), x, y, z);

end

