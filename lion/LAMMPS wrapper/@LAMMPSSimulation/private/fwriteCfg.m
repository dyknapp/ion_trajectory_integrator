function fwriteCfg( fileHandle, s )
%FWRITECFG Writes a struct of cfg file handles to the given file.

if ~isfield(s, 'cfgFileHandle')
    error('Invalid input');
end

for i=1:length(s)
    fwriteCell(fileHandle, s(i).cfgFileHandle());
end

end
