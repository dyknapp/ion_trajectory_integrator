function [ strings ] = langevinBathCfg( fixID, temperature, dampingTime, group, seed )
%LANGEVINBATHCFG Creates the cfg file txt for a langevin bath

if nargin < 5
    seed = 1337;
end

if nargin < 4
    groupName = 'all';
else
    if ~isa(group, 'LAMMPSGroup')
       error('group variable must be of type LAMMPSGroup'); 
    end
    
    groupName = group.ID;
end

strings = { '#Adding a langevin bath...' };
strings{end+1} = sprintf('fix %s %s langevin %e %e %e %d', fixID, groupName, temperature, temperature, dampingTime, seed); 
strings{end+1} =  '';

end

