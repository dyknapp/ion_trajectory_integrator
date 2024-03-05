function Unfix( sim, fix )
% UNFIX Disables a previously added fix. 
%
% Unfix inserts a command to the simulation input file that disables a fix
% previously added to the simulation via AddFix.
%
% SYNTAX: Unfix(sim, fix)
%
% Example:
%  Unfix(sim, biasField)
% See Also: Add, Remove

% validate input
if ~isa(fix, 'LAMMPSFix')
    error('Input argument invalid; must be of class LAMMPSFix.');
end

unfix = InputFileElement();
unfix.createInputFileText = @(~) unfixCfg(getID(fix));
sim.Add(unfix);

end

