function [ id ] = getUnusedID( type )
%GETUNUSEDID Gets an unused ID. Valid types are 'fix', 'region',
%'atomtype', 'group', 'priority', 'command'. All return integer except fix, which returns a
%string.
% This command is for internal use of the lammps wrapper, in particular for
% generating unique ids of fixes, groups etc when writing the input file
% for LAMMPS. This should not be used by the end user.

persistent currAtomTypeID currFixID currRegionID currGroupID currPriority ...
    currCommID currDumpID;

if ~exist('currAtomTypeID', 'var') || size(currAtomTypeID,1) == 0
    currAtomTypeID = 1;
    currFixID = 1;
    currRegionID = 1;
    currGroupID = 1;
    currPriority = 1;
    currCommID = 1;
    currDumpID = 1;
end

switch type
    case 'fix' %requires a string identifier, which we generate in a deterministic way
        id = currFixID;
        currFixID = currFixID + 1;
        
        %convert to a string we can use. We count in terms of aa->az,
        %ba->bz etc
        id = getStringID(id);
    case 'region'
        id = currRegionID;
        currRegionID = currRegionID + 1;
    case 'atomtype'
        id = currAtomTypeID;
        currAtomTypeID = currAtomTypeID + 1;
    case 'group'
        id = currGroupID;
        currGroupID = currGroupID + 1;
    case 'reset'
        clear currAtomTypeID currFixID currRegionID currGroupID;
        id = 0;
    case 'priority'
        id = currPriority; %used for ordering elements inside lammps input
        currPriority = currPriority + 1;
    case 'command'
        id = currCommID;
        currCommID = currCommID + 1;
    case 'dump'
        id = getStringID(currDumpID);
        currDumpID = currDumpID + 1;
end
end

