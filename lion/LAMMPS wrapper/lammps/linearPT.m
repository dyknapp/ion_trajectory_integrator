function [ fix ] = linearPT( OscillatingV, EndcapV, Z0, R0, GeometricConstant, RFFrequency, Anisotropy, offset )
%LINEARPT Applies an oscillating electric field to atoms. The
%characterisation of the trap follows the parameters outlined in Berkeland
%et al. (1998) for a linear Paul trap geometry.
% RFFrequency is specified in Hz. Z0 and R0 are the axial and radial
% dimensions of the trap. The optional parameter Anisotropy is used to
% imbalance fields in X and Y directions, such that V_y = anisotropy * V_x
%
% This method evaluates the force due to the full radiofrequency electric
% fields - pseudopotential approximations may be more appropriate for most
% cases.
%
% OscillatingV and RFFrequency may be specified as vectors, in which case a
% multi-frequency trap is created as discussed in
% https://arxiv.org/pdf/1310.6294.pdf
%
% SYNTAX: linearPaulTrap(OscV, endV, z0, r0, geomC, RFFreq, Anisotropy, offset)
% 
% See Also: LinearSecularPaulTrap, http://tf.nist.gov/general/pdf/1226.pdf

if nargin < 8
    offset = [0 0]';
elseif length(offset) ~= 2
    error('Offset must a 2d vector');
end

if nargin < 7
    Anisotropy = 1;
end

if length(OscillatingV) ~= max(size(OscillatingV))
    error('OscillatingV must be a vector or scalar.');
end

if length(RFFrequency) ~= max(size(RFFrequency))
    error('RFFrequency must be a vector or scalar.');
end

RFFrequency = RFFrequency(:);
OscillatingV = OscillatingV(:);

fix = LAMMPSFix();
fix.time = 1/max(RFFrequency(:))/20;
fix.createInputFileText = @linearPaulTrapCfg;
fix.InputFileArgs = { OscillatingV, EndcapV, Z0, R0, GeometricConstant, RFFrequency, Anisotropy, offset };

end