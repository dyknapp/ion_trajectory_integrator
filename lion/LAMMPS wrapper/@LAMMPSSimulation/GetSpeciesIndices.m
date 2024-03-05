function indices = GetSpeciesIndices(sim)
%GETSPECIESINDICES Returns a cell array. Each element of the array is a
%list of atomic indices of a particular atom type. The array is ordered in
%the same order that atom types are defined using AddAtomType.
%
% SYNTAX: GetSpeciesIndices(sim)
%
% Example:
%  i = GetSpeciesIndices(sim)
%  calcium = i{1};
% See Also: AddAtomType, AddAtoms

indices = cell(length(sim.AtomTypes), 1);
atoms = sim.AtomList;
for i=1:length(sim.AtomTypes)
    type = sim.AtomTypes(i);
    atomTypes = [atoms.Type];
    mask = [atomTypes.ID] == type.ID;
    subset = atoms(mask);
    indices{i} = [subset.ID];
end

end