function id = AddAtom( sim, atom )
%ADDATOM Adds an atom to the ion trap simulation.
% The atom should be of type AtomPlacement. This function is automatically
% called when an instance of AtomPlacement is created. To add large numbers
% of atoms to a simulation, see createIonCloud or placeAtoms.
% 
% SYNTAX: AddAtom(sim, atom)
% 
% See Also: CREATEIONCLOUD, PLACEATOMS

sim.assertNotRun();

if ~isa(atom, 'AtomPlacement')
    error('atom must be of type AtomPlacement');
end

if numel(atom) ~= 1
    error('Atoms must be added one at a time.');
end

id = length(sim.AtomList)+1;
sim.AtomList(id) = atom;

end
