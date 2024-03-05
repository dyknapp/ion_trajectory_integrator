function [ strings ] = linearPaulTrapCfg(fixID, OscillatingV, EndcapV, Z0, R0, GeometricConstant, RFFrequency, Anisotropy, offset )
%LINEARPAULTRAPCFG Writes the config file for a linear paul trap.

% Displacement of the trap centre:
xPos = sprintf('(x-%e)', offset(1));
yPos = sprintf('(y-%e)', offset(2));
if offset(1) == 0; xPos = 'x'; end
if offset(2) == 0; yPos = 'y'; end

blockTitle = sprintf('#Creating a Linear Paul Trap... (fixID=%s)', fixID);

% Calculate potential gradients
Vgrad = OscillatingV / (R0).^2;
Vstat = GeometricConstant * EndcapV / (Z0).^2;

% Create strings describing voltages.
voltStrx = {};
voltStry = {};
for i=1:length(Vgrad)
    voltStrx{i} = sprintf('%s * cos(%s * step*dt)', Vgrad(i), RFFrequency(i)*2*pi );
    voltStry{i} = sprintf('%s * cos(%s * step*dt)', Vgrad(i)*Anisotropy, RFFrequency(i)*2*pi );
end
voltStrx = strjoin(voltStrx, '+');
voltStry = strjoin(voltStry, '+');

% Calculate per-atom forces
Ex_varname = [ 'Ex' fixID ];
Ey_varname = [ 'Ey' fixID ];
Ez_varname = [ 'Ez' fixID ];
ExStr = sprintf('variable %s atom "(%s) * (%s) + %.8f * x"', Ex_varname, voltStrx, xPos, Vstat);
EyStr = sprintf('variable %s atom "(%s) * -(%s) + %.8f * y"', Ey_varname, voltStry, yPos, Vstat);
EzStr = sprintf('variable %s atom "-%.8f * z"', Ez_varname, 2*Vstat);

% Add fix
fixStr = sprintf('fix %s all efield v_%s v_%s v_%s', fixID, Ex_varname, Ey_varname, Ez_varname);
strings = { blockTitle, ExStr, EyStr, EzStr, fixStr };

end

