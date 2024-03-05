function [ strings ] = nveIntegratorCfg( fixID, groupName, preamble )
%NVEINTEGRATORCFG Configures all atoms to undergo nve integration.

if nargin < 2
    groupName = 'all';
end

strings = { '#Configuring integration for all atoms', preamble, sprintf('fix %s %s nve', fixID, groupName), ''};

end

