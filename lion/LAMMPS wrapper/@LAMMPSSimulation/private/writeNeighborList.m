function writeNeighborList( sim, fHandle )
%WRITENEIGHBORLIST Writes text to the config file that configures the
%neighbor list calculation.
% See: http://lammps.sandia.gov/doc/neigh_modify.html
%      http://lammps.sandia.gov/doc/neighbor.html

fprintf(fHandle, '# Configure Neighbor list: \n');
fprintf(fHandle, 'neighbor %.3f %s\n', sim.NeighborSkin, sim.NeighborList);
fprintf(fHandle, 'neigh_modify once yes\n\n');

end

