function writeGroups( sim, fHandle )
%WRITEGROUPS Writes group definitions to the output file

fprintf(fHandle, '# Group definitions:\n');
fprintf(fHandle, sprintf('# (%d group(s))\n', length(sim.Groups)));

for i=1:length(sim.Groups)
    g = sim.Groups(i);
    fprintf(fHandle, g.getInputFileText());
end
fprintf(fHandle, '\n');


end

