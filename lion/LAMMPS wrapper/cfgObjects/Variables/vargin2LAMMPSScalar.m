function [ vars ] = vargin2LAMMPSScalar( cellin )
%vargin2LAMMPSScalar Returns a cell for the given input cell that
%represents a list of symbols that represent scalars in LAMMPS.

if ~iscell(cellin)
    c = cell(1);
    c{1} = cellin;
    cellin = c;
end

%Examine form of variables: is it a fix/variable we need to look at or is
%it a 'special case' of inbuilt lammps variables, eg x, y, z
vars = cell(0);

for j=1:length(cellin)
    variable = cellin{j};
    
    if ~isa(variable,'LAMMPSVariable')
        
        if ~isValidLAMMPSVariable(variable)
            error('Not a valid LAMMPS Variable: %s', variable);
        end
        
        vars{end+1} = variable;
        
    else
        %concatenate outputs from variable
        vars = [vars variable.getOutputs()];
    end
end

end

