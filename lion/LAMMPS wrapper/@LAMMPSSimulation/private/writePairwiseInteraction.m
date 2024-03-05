function writePairwiseInteraction( sim, fHandle)
%WRITEPAIRWISEINTERACTION Configures pairwise interactions between atoms.
%Sets scaling of Coulomb repulsion to 1 for all atom types.

fprintf(fHandle, '# Configure pairwise interactions for long range Coulombics only\n');
fprintf(fHandle, 'pair_style coul/cut %e\npair_coeff * * \n\n', sim.CoulombCutoff);


end

