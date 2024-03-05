function [ out ] = lionVersion( )
%LIONVERSION Gets the current version of the lion wrapper

curr = cd;
out = '???';

a = mfilename('fullpath');
cd(fileparts(a));
try
    [~,str] = system('git rev-parse HEAD');
    out = strtrim(str);
catch
    warning('Could not determine version of repository. Is git installed and on system path?');
end
cd(curr);

end