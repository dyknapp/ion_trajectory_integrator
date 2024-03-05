function [ strings ] = minimizeCfg( etol, ftol, maxiter, maxeval, maxDist )
%MINIMIZECFG Minimise energy via atom positions
strings =           { '#minimize', 'min_style quickmin', sprintf('min_modify dmax %e', maxDist), sprintf('minimize %e %e %d %d', etol, ftol, maxiter, maxeval )};
end

