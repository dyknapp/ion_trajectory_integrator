function [ outLV ] = stepAvg( lammpsVariables, nSteps )
%STEPAVG Averages the given cell of lammpsVariables over nSteps. nSteps may
%be a zero-argument function handle that returns an integer > 0, in which case evaluation is delayed
%until the input file is written.
% Output: a 1x1 cell of LAMMPSVariable
% Example:
%  dump('output.txt', {'id',stepAvg({'x','y','z'}, 10)}, 1);
% See Also: timeAvg

%Potential casue of error: stepAvg does not respect when LAMMPSVariables
%are defined, eg step averaging a quantity that is defined only every n
%frames will throw an error as it attempts to query the value of this every
%frame.

if isempty(lammpsVariables)
    error('lammpsVariables cannot be empty');
end

%input can be either a single literal, lammps variable, or a cell of either
%of these.
if ~iscell(lammpsVariables) && ~(isa(lammpsVariables, 'LAMMPSVariable') || ischar(lammpsVariables))
    error('lammpsVariables input must either be a string literal, a LAMMPSVariable object or a (mixed) cell array of these.');
end

%convert single property to length 1 cell.
if ~iscell(lammpsVariables)
    lv = cell(1,1);
    lv{1} = lammpsVariables;
    lammpsVariables = lv;
end

%Convert atomic property literals to lammps properties.
lammpsVariables = literalsToLammpsVariables(lammpsVariables);

%Cannot time average mass, id, charge, atom type...
for lv=lammpsVariables
    if ~canTimeAverage(lv{1})
        error('Quantity is not suitable for time averaging: %s', lv{1}.Text);
    end
end

outLV = LAMMPSVariable();
outLV.Frequency = nSteps;
outLV.setVariableType('fix');
outLV.Length = sum(cellfun(@(x) x.Length, lammpsVariables));

%recursively inherit preambles of input lammps variables eg to define an
%other needed quantities.
preamble = {};
for lv=lammpsVariables
    strs = getLines(lv{1});
    if ~isempty(strs)
        for s=strs
            preamble{end+1} = s;
        end
    end
end

%concatenate outputs for each variable
variableString = '';
for lv=lammpsVariables
    variableString = sprintf('%s %s', variableString, strjoin(lv{1}.getOutputs()));
end

%Add preamble specific to time averaging. Note that we have to delay
%evaluation if nSteps is a function handle until the input file is written,
%eg to time average over correct number of steps for duration.

%convert to function handle
if isa(nSteps, 'function_handle')
    nStepsHandle = nSteps;
else
    nStepsHandle = @() nSteps;
end

preamble{end+1} = '# define time integration of variables over number of steps:';
outLV.createInputFileText = @(~) [preamble, ...
    {sprintf('fix %s all ave/atom 1 %d %d%s', getID(outLV), nStepsHandle(), nStepsHandle(), variableString)}];

end

