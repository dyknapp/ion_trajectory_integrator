function Add( sim, inputFileElement )
% ADD Adds a fix or command to the lammps simulation. The input argument is
% the fix or command to add to the simulation.
%
% SYNTAX: Add(obj)
%
% Example:
%  Add(InputFileElement)
%
% See Also: DUMP, MINIMIZE, RUNCOMMAND, THERMALVELOCITIES,
% LINEARPT, SHO, LANGEVINBATH, LASERCOOL, EFIELD, HARMONICOSCILLATOR

sim.assertNotRun();

if isa(inputFileElement, 'InputFileElement')
    sim.Elements{end+1} = inputFileElement;
    
    % Assign the input element a unique id
    inputFileElement.ID = idString(6);
    
    if isa(inputFileElement, 'LAMMPSFix')
        fix = inputFileElement;
        if fix.time > 0 && fix.time < sim.LimitingTimestep
            sim.LimitingTimestep = fix.time;
        end
    end
else
    error('Invalid input argument: Input argument must be either a LAMMPSFix or valid subclass of PrioritisedCfgObject (eg LAMMPSRunCommand, LAMMPSDump)');
end

end