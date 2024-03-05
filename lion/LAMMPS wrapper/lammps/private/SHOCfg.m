function [ strings ] = SHOCfg(fixID, k_x, k_y, k_z, groupID )
%SHOCFG Config file for atoms attached to the origin by springs. GroupID
%should be a string, or leave blank for 'all'

strings =           {'#SHO'};

if (nargin < 4)
    groupID = 'all'; %Fix for all atoms if no atomType is specified.
end

strings{end+1} =     sprintf('variable k_x%s\t\tequal %e', fixID, k_x);
strings{end+1} =     sprintf('variable k_y%s\t\tequal %e', fixID, k_y);
strings{end+1} =     sprintf('variable k_z%s\t\tequal %e', fixID, k_z);
strings{end+1} =     sprintf('variable fX%s atom "-v_k_x%s * x"', fixID, fixID);
strings{end+1} =     sprintf('variable fY%s atom "-v_k_y%s * y"', fixID, fixID);
strings{end+1} =     sprintf('variable fZ%s atom "-v_k_z%s * z"', fixID, fixID);
%Set up energy for use in minimisation
strings{end+1} =	 sprintf('variable E%s atom "v_k_x%s * x * x / 2 + v_k_y%s * y * y / 2 + v_k_z%s * z * z / 2"', fixID, fixID, fixID, fixID);
strings{end+1} =     sprintf('fix %s %s addforce v_fX%s v_fY%s v_fZ%s energy v_E%s', fixID, groupID, fixID, fixID, fixID, fixID);
strings{end+1} = '';
end

