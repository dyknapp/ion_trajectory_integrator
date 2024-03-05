function [ fix ] = linearPseudoPT( OscillatingV, EndcapV, Z0, R0, GeometricConstant, RFFrequency, atomType )
%LINEARPSEUDOPT Pseudopotential approximation of a linear paul trap.
% This approximation replaces the radiofrequency electric field of the
% linear Paul trap with a harmonic pseudopotential. This harmonic
% pseudopotential may be approximated via the Mathieu equation (see eg.
% Berkeland's paper, J. App. Phys. 83, (10) 1998).
%
% As the pseudopotential is dependent on the charge:mass ratio of the ion,
% this fix requires that an atomType be supplied. If multiple atom types
% exist in the trap, one fix should be created for each type.
% 
% SYNTAX: linearPseudoPT( OscV, endcapV, z0, r0, geomC, RF, atomType )
% 
% See Also: linearPaulTrap, http://tf.nist.gov/general/pdf/1226.pdf

a_r = -4 * (atomType.Charge * 1.6e-19) * GeometricConstant * EndcapV / ...
    ((atomType.Mass *  1.660e-27) * Z0^2 * (RFFrequency * 2 * pi)^2);
a_z = - 2 * a_r;

q_r = 2 * (atomType.Charge * 1.6e-19) * OscillatingV / ...
    ((atomType.Mass *  1.660e-27) * R0^2 * (RFFrequency * 2 * pi)^2);

f_r = ((RFFrequency * 2 * pi) / 2 * (a_r + (q_r^2)/2)^0.5) / 2 / pi;
f_z = ((RFFrequency * 2 * pi) / 2 * (a_z)^0.5 ) / 2 / pi;

fprintf('Frequency of motion: f_r=%e, f_z=%e\n', f_r, f_z)

fix = harmonicOscillator(f_r, f_r, f_z, atomType.Group);

end

