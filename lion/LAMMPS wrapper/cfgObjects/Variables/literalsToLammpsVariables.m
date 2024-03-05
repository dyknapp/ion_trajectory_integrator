function [lammpsVariables] =  literalsToLammpsVariables(lammpsVariables)
%literalsToLammpsVariables Converts any literal strings in a mixed cell of
%strings and LAMMPSVariables into LAMMPSVariables. Returns a pure cell of
%LAMMPSVariables.

for i=1:length(lammpsVariables)
    
    if ischar(lammpsVariables{i})
        switch lammpsVariables{i}
            case {'x','y','z','vx','vy','vz','mass','q','id','type'}
                
                %Create a lammpsvariable object to replace literal in
                %cell array.
                lv = LAMMPSVariable;
                lv.Frequency = @() 1; %literals always defined.
                lv.setVariableType('atomProperty');
                lv.Text = lammpsVariables{i};
                lammpsVariables{i} = lv;
                
            otherwise
                error('Invalid string literal for lammps variable. Valid properties are x, y, z, vx, vy, vz, mass, q, id, type.');
        end
    end
    
    %literals are now replaced.
    if ~isa(lammpsVariables{i}, 'LAMMPSVariable')
        error('lammpsVariables input must either be a string literal, a LAMMPSVariable object or a (mixed) cell array of these.');
    end
end
end