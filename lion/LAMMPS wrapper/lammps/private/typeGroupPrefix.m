function [ strings ] = typeGroupPrefix( cfgFileHandleWithGrpArg, LAMMPSatomType )
%TYPEGROUPPREFIX Prefixes the supplied CFG text with lammps input to create
%a group from an atom type.
% The supplied cfg file handle must take group as an argument.

%Default to all.
groupID = 'all';

if exist('LAMMPSatomType', 'var')
    strings = { '#Prefix: generate group for atom type.' };
    
    if (nargin < 3)
        groupID = sprintf('group%d', getUnusedID('group'));
        %create a group for this atom type we can use to apply the fix.
        strings{end+1} = sprintf('group %s type %d', groupID, LAMMPSatomType);
    end
else
    strings = { ' ' };
end

strings = [strings cfgFileHandleWithGrpArg(groupID)];

end