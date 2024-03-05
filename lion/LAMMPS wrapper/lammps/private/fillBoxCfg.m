function [ strings ] = fillBoxCfg( radius, charge, amuMass, number, spacing )
%FILLBOXCFG Writes the config file to fill a described box with atoms.

%get required ids.
regionID = getUnusedID('region');
atomTypeID = getUnusedID('atomtype');

strings =           {'# Creating Ion Cloud....'};
strings{end+1} =    sprintf('lattice fcc %f', spacing);
strings{end+1} =    sprintf('region %d sphere 0.0 0.0 0.0 %.8f', regionID, radius);
strings{end+1} =    sprintf('create_atoms %d random %d 1337 %d', atomTypeID, number, regionID);
strings{end+1} =    sprintf('mass %d %e', atomTypeID, 1.660e-27 * amuMass);
strings{end+1} =    sprintf('set type %d charge %e', atomTypeID, 1.6e-19 * charge);
strings{end+1} =    '';
end

