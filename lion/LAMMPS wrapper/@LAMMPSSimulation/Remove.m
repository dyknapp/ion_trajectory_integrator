function Remove( sim, obj )
% REMOVE Disables a previously added fix or dump (obj) in the lammps
% simulation. This can be used to eg. suspend output to a file or turn off
% a force previously applied to the atoms.
%
% SYNTAX: Remove(sim, obj)
%
% Example:
%  Remove(LAMMPSFix)
%  Remove(LAMMPSRunCommand)
% See Also: Add, Unfix

if isa(obj, 'LAMMPSFix')
    sim.Unfix(obj);
elseif isa(obj, 'LAMMPSDump')
    rc = InputFileElement();
    rc.createInputFileText = @ (~) { '#Remove dump', sprintf('undump %s', obj.getID())};
    sim.Add(rc);
else
    error('Invalid argument to LAMMPSSimulation.Remove: Must be either a LAMMPSFix or a LAMMPSDump');
end

end

