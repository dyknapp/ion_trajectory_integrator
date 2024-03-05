function atomType = AddAtomType( sim, charge, amuMass )
%ADDATOMTYPE Define atomic species used in simulations.
% Units are in atomic charges and mass units. Returns a reference to the
% created atom type, which can be used with atom placement commands.
%
% SYNTAX: AddAtomType(sim, charge, amuMass)
%
% Example: 
%  AddAtomType(sim, 1, 40) % calcium 40+ ion
% See Also: AddAtoms

sim.assertNotRun();

id = length(sim.AtomTypes(:))+1;
atomType = AtomType(sim, id, charge, amuMass);
sim.AtomTypes(end+1) = atomType;

end
