function [ strings ] = harmonicOscillatorCfg(fixID, f_x, f_y, f_z, group )
%SHOCFG Creates a harmonic oscillator centred around the origin.
% If a group is specified the fix will apply to only those atoms.

if (nargin < 4)
    groupName = 'all'; %Fix for all atoms if no atomType is specified.
else
    if ~isa(group, 'LAMMPSGroup')
       error('group must be a LAMMPSGroup'); 
    end
    
    groupName = group.ID;
end

strings =           {'# Harmonic Oscillator: '};
strings{end+1} =     sprintf('variable %s_w_x\t\tequal %e', fixID, (2*pi*f_x).^2);
strings{end+1} =     sprintf('variable %s_w_y\t\tequal %e', fixID, (2*pi*f_y).^2);
strings{end+1} =     sprintf('variable %s_w_z\t\tequal %e', fixID, (2*pi*f_z).^2);
strings{end+1} =     sprintf('variable %s_fX atom "-v_%s_w_x * x * mass"', fixID, fixID);
strings{end+1} =     sprintf('variable %s_fY atom "-v_%s_w_y * y * mass"', fixID, fixID);
strings{end+1} =     sprintf('variable %s_fZ atom "-v_%s_w_z * z * mass"', fixID, fixID);

%Set up energy for use in minimisation:
strings{end+1} =	 sprintf('variable %s_E atom "v_%s_w_x * x * x / 2 + v_%s_w_y * y * y / 2 + v_%s_w_z * z * z / 2"', fixID, fixID, fixID, fixID);

strings{end+1} =     sprintf('fix %s %s addforce v_%s_fX v_%s_fY v_%s_fZ energy v_%s_E', fixID, groupName, fixID, fixID, fixID, fixID);
strings{end+1} = '';
end

