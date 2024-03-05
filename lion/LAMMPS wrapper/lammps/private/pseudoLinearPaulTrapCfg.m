function [ strings ] = pseudoLinearPaulTrapCfg(fixID, k_r, k_z, atomType )
%LINEARPAULTRAPCFG Config file for atoms attached to the origin by springs,
%cylindrical symmetry.

%TODO: Change this to use type group prefixing.

strings =           {sprintf('#Pseudopotenial approximation for Linear Paul trap... (fixID=%s)', fixID)};

if (nargin < 3)
   error('Pseudopotential approximation must be applied individually to each atom type.');
else
   groupID = sprintf('group%d', getUnusedID('group'));
   %create a group for this atom type we can use to apply the fix.
   strings{end+1} = sprintf('group %s type %d', groupID, atomType);
end

%Add a cylindrical SHO for the pseudopotential.
temp = SHOCfg(fixID, k_r, k_r, k_z, groupID);
for i=1:length(temp)
   strings{end+1} = temp{i}; 
end

end

