function [ n ] = gasDensity( T, mbar)
%GASDENSITY Calculate the number density of particles in an ideal gas.

%1mbar = 100Pa
kb = 1.3806488e-23;
n = mbar*100/(T*kb);
end