function [ strings ] = thermalVelocitiesCfg( temperature, zeroTotalMom, seed )
%THERMALVELOCITIESCFG generates config file text to set all atoms to a
%given temperature.

strings = { sprintf('velocity all create %e %d mom %s rot yes dist gaussian', temperature, seed, zeroTotalMom)};

end

