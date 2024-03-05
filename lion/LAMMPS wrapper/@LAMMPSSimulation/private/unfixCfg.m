function [ strings ] = unfixCfg( fixID )
%DELETEFIXCFG Config file text to remove a fix.

strings = {'#Deleting a fix', sprintf('unfix %s', fixID)};
end

