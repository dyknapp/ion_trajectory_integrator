function [ strings ] = laserCoolCfg(fixID, groupID, gamma )
%LASERCOOLCFG Config file excerpt for laser cooling of an atom species.
% Gamma encodes strength of cooling in norm(gamma) and the direction of
% cooling in gamma/norm(gamma)

coolDir = gamma / norm(gamma);
coolForce = norm(gamma);

strings = {'#Define laser cooling for a particular atom species.'};

% Define variable corresponding to $ -v \cdot coolDir $
strings{end+1} =     sprintf('variable %s_vel atom "%e * vx + %e * vy + %e * vz"', fixID, coolDir(1), coolDir(2), coolDir(3));

% strings{end+1} =     sprintf('variable %s_gammaX\t\tequal %e', fixID, coolDir(1).^2 * coolForce);
% strings{end+1} =     sprintf('variable %s_gammaY\t\tequal %e', fixID, coolDir(2).^2 * coolForce);
% strings{end+1} =     sprintf('variable %s_gammaZ\t\tequal %e', fixID, coolDir(3).^2 * coolForce);

strings{end+1} =     sprintf('variable %s_fX atom "-v_%s_vel * mass * %e"', fixID, fixID, coolForce * coolDir(1));
strings{end+1} =     sprintf('variable %s_fY atom "-v_%s_vel * mass * %e"', fixID, fixID, coolForce * coolDir(2));
strings{end+1} =     sprintf('variable %s_fZ atom "-v_%s_vel * mass * %e"', fixID, fixID, coolForce * coolDir(3));
strings{end+1} =     sprintf('fix %s %s addforce v_%s_fX v_%s_fY v_%s_fZ', fixID, groupID, fixID, fixID, fixID);

strings{end+1} =     '';

end