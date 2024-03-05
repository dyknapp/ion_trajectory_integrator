function [ outLV ] = timeAvg( lammpsVariables, time )
%TIMEAVG Periodically averages a cell of LAMMPSVariable (or literals representing
%LAMMPSVariable) over the specified duration.
% Output: a 1x1 cell of LAMMPSVariable
% Example:
%  dump('secularV.txt', {'id',timeAvg({'vx','vy','vz'}, 1/RF)}, 1);
% See Also: stepAvg

%Create function handle to calculate correct number of timesteps for
%required duration. This will be invoked when writing the input file.

if (time <= 0)
    error('time must be a positive quantity greater than zero.');
end

nStepsHandle = @() ceil(time/cfghelperTimestep('get'));
outLV = stepAvg(lammpsVariables, nStepsHandle);

end

