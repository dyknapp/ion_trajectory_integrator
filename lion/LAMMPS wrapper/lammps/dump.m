function [ dumpCmd ] = dump( filename, variables, interval )
%DUMP Dumps variables from lammps into files for analysis.
% Outputs the given variables into the specified file every internal
% simulation steps. variables is a cell of LAMMPSVariables and/or
% the literals x,y,z,vx,vy,vz. The interval defaults to 10 if not specified.
%
% SYNTAX: dump( filename, variables, interval )
%
% Example:
%  dump('output.txt', {'id','x','y','z'}, 10);
% See Also: LAMMPSAverage

if nargin < 3
    interval = 10;
end

if isempty(variables)
    error('lammpsVariables cannot be empty');
end

%input can be either a single literal, lammps variable, or a cell of either
%of these.
if ~iscell(variables) && ~(isa(variables, 'LAMMPSVariable') || ischar(variables))
    error('lammpsVariables input must either be a string literal, a LAMMPSVariable object or a (mixed) cell array of these.');
end

%convert single property to length 1 cell.
if ~iscell(variables)
    lv = cell(1,1);
    lv{1} = variables;
    variables = lv;
end

%Convert atomic property literals to lammps properties.
variables = literalsToLammpsVariables(variables);

%If no error thrown, lammpsVariables is now a cell of LAMMPSVariable objects.
dumpCmd = LAMMPSDump();
dumpCmd.createInputFileText = @dumpCfg;
dumpCmd.InputFileArgs = { filename, variables, interval };

end
