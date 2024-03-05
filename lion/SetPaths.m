if isempty(strfind(path,'LAMMPS Wrapper'))
    addpath(genpath(sprintf('%s\\LAMMPS Wrapper', cd)))
    addpath(genpath(sprintf('%s\\Analysis', cd)))
    addpath(genpath(sprintf('%s\\DCDplotting', cd)))
    addpath(genpath(sprintf('%s\\Tools', cd)))
    addpath(genpath(sprintf('%s\\Util', cd)))
    addpath(genpath(sprintf('%s\\Plot', cd)))
end