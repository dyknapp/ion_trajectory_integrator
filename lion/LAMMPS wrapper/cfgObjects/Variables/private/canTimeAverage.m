function [ result ] = canTimeAverage( lammpsVariable )
%CANTIMEAVERAGE Attempts to estimate whether the given lammpsVariable is
%valid for time averaging. Can only assert the property can't be time
%averaged (false), not that it can.

%Note: this probably won't work for second-order variables, eg time
%averages of other time averages. However, I include the feature so that in
%the likely cases to cause issue (eg attempting to time average an atom id)
%it will confront the user.

if ~isa(lammpsVariable, 'LAMMPSVariable')
    error('lammpsVariable must be an object of class LAMMPSVariable');
end

if strcmp(lammpsVariable.variableType, 'atomProperty')
    switch lammpsVariable.Text
        case {'id', 'mass', 'charge', 'type'}
            result = false;
        otherwise
            result = true;
    end        
else
    result = true; %can't tell for this lv.
end

end

