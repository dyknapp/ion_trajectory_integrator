function [ s ] = getVariablePrefix( type )
%GETVARIABLETYPE Gets the variable prefix suitable for the given variable
%type.
switch type
    case 'fix'
        s = 'f_';
        
    case 'var'
        s = 'v_';
        
end

end

