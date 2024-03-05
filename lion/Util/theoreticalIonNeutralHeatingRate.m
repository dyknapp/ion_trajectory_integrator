function [ h_coll, gamma_elastic ] = theoreticalIonNeutralHeatingRate( n, Tn, mi, p, Ti, mn)
%THEORETICALHEATINGRATE Calculates theoretical heating rate for an ion
%cloud due to ion-neutral collisions.
% This is based on eqn 8 from Zhang et. al Phys Rev A 2007
% n - number density of neutral particles.
% Tn - temperature of neutral background
% mi - mass of ion
% Optional parameters:
% p - polarizability of neutral particle (default: 1.76e-34, value for N2)
% Ti - temperature of ion cloud (default: 0);
% mn - mass of neutral particles

if nargin < 3
   error('Not enough input parameters'); 
end

if ~exist('p', 'var')
   p = 1.76e-34*1e-6; 
end

if ~exist('mn', 'var')
   mn = 2.32587e-26*2;
end

if ~exist('Ti', 'var')
   Ti = 0;
end

%reduced mass
u = (mn*mi)/(mi+mn);

h_coll = (3*2.21/4)*(1.6021766e-19)/(8.854187e-12)*(1.3806488e-23)*(u.^0.5/(mi+mn))*(Tn-Ti)*(n*p.^0.5);
 gamma_elastic =2.21/4*(1.6e-19)/(8.854187e-12)*n*(p/u).^0.5;
end