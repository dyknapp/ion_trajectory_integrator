function [ oscV, endV ] = getVs4aq( ionType, RFFrequency, z0, r0, geometricC, a, q )
%GETVS4AQ Calculates the oscillating and endcap voltages required for
%specified Mathieu (a,q) stability parameters in the defined trap.
% Assumes a linear Paul trap geometry as defined in Berkeland et. al.
% (1998). ionType is the species to calculate voltages for. a and q
% parameters are specified for the x-axis motion.

ionMass = ionType.Mass * 1.66e-27;
ionCharge = ionType.Charge * 1.6e-19;

endV = a * ionMass * z0^2 * (2*pi*RFFrequency)^2/-(geometricC*4*ionCharge);
oscV = q * -(ionMass * r0^2 * (2*pi*RFFrequency)^2) / (2*ionCharge);

end