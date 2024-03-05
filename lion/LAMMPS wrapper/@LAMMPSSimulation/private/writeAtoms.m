function writeAtoms( sim, fHandle )
%WRITEATOMS Writes information describing the initial atom ensemble to the
%log file

%Write atom insertion commands
fprintf(fHandle, '# Placing atoms\n');
for i=1:length(sim.AtomList)
    atom = sim.AtomList(i);
    fprintf(fHandle, 'create_atoms %d single %e %e %e units box\n', atom.Type.ID, atom.X, atom.Y, atom.Z);
end

%Definitions for atom types
fprintf(fHandle, '# Atom definitions\n');
for i=1:length(sim.AtomTypes)
    atomType = sim.AtomTypes(i);
    fprintf(fHandle, 'mass %d %e\n', atomType.ID, 1.660e-27 * atomType.Mass);
    fprintf(fHandle, 'set type %d charge %e\n', atomType.ID, 1.6e-19 * atomType.Charge);
end
fprintf(fHandle, '\n');

end