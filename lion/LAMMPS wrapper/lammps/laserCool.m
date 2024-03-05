function [ fixObj ] = laserCool( atomType, gamma )
%LASERCOOL A simplistic model of laser cooling.
% Simulates the laser cooling of atoms along a given axis by applying an
% anisotropic viscous damping force. The laser cooling is applied only to
% the specified atomType. gamma = [gx, gy, gz] defines the strength of the
% damping force of the form f_i = - gamma_i m * v_i. gamma is equivalent to
% 1/damping time.
% 
% Syntax: laserCool( atomType, gamma )
%
% Example: laserCool( calcium, [ 1 Inf Inf ] / 60e-6 ) % 60us viscous force along one direction
% 
% See Also: langevinBath

if ~isa(atomType, 'AtomType')
   error('must specify a valid atom species.'); 
end

if length(gamma) == 1
    repmat(gamma,1,3);
elseif length(gamma) ~= 3
   error('gamma must be a length-3 vector or a single number'); 
end

if any(gamma < 0)
   error('damping must be greater than zero'); 
end

if any(isinf(gamma) | isnan(gamma))
    error('invalid damping parameter: must be finite and numerical');
end

fixObj = LAMMPSFix();
fixObj.time = 0.1 / max(gamma(:));
fixObj.createInputFileText = @laserCoolCfg;
fixObj.InputFileArgs =  { atomType.Group.ID, gamma };

end

