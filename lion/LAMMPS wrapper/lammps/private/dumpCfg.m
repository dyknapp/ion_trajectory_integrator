function [ strings ] = dumpCfg( id, filename, lammpsVariables, Nsteps )
%DUMPCFG Writes cfg file text for dumping lammps variables to file
%lammpsVariables is a cell of LAMMPSVariable objects.

strings = { '#Dump output' };

%determine base output frequency of data as LCM of all lv.frequency. This
%ensures all variables are defined when the dump is written.
baseFreq = 1;
for lv=lammpsVariables
   baseFreq = lcm(baseFreq, lv{1}.Frequency()); 
end

outputFrequency = baseFreq * Nsteps;

vars = vargin2LAMMPSScalar(lammpsVariables);
outputVars = sprintf('%s ', vars{:});

%Write preamble for lammpsVariables - some lammpsVariables may need to
%define some other explicit text, eg a lammpsVariable produced by timeAvg
%may need to setup the initial time averaging fix. Individual
%lammpsVariable styles should implement recursion if it is needed (eg
%stepAvg of a timeAvg)
for i=1:length(lammpsVariables)
   extra = getLines(lammpsVariables{i});
   if ~isempty(extra)
       for str=extra
          strings{end+1} = str{1};
       end
   end
end

strings{end+1} =  sprintf('dump %s all custom %d %s %s', id, outputFrequency, filename, outputVars);

end

