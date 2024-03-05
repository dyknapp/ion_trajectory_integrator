function [ strings ] = efieldCfg( fixID, Ex, Ey, Ez )
%EFIELDCFG Writes cfg file text for a constant efield.

strings = { '#Static E-field', sprintf('fix %s all efield %e %e %e', fixID, Ex, Ey, Ez) };

end

