function writeTimeIntegration(sim, fHandle)

fprintf(fHandle, '# Time integration: \n');
% If 'Rigid' fixes exist we need to separate them from the typical nve
% integrator.
nveGroup = 'group nonRigidBody union all';

if any([sim.Groups.Rigid])
    rigidGroups = find([sim.Groups.Rigid]);
    
    % construct rigid string representing rigid bodies
    rigidStr = [];
    for i=1:length(rigidGroups)
        rigidStr = [ rigidStr sim.Groups(rigidGroups(i)).ID ' '];
    end
    
    nveGroup = ['group nonRigidBody subtract all ' rigidStr];
    fprintf(fHandle, ['group rigidAtoms union ' rigidStr '\n']);
    fprintf(fHandle, 'fix rigidIntegrator rigidAtoms rigid group %d %s\n', length(rigidGroups), rigidStr);
end

fprintf(fHandle, '%s\n', nveGroup);
fprintf(fHandle, 'fix timeIntegrator nonRigidBody nve\n');

end