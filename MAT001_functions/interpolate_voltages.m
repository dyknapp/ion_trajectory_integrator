function interpolated_voltages = interpolate_voltages(voltages, t, step_times)
% A function finding the electrode voltage at an arbitrary time by
% interpolating based on the input parameter voltages, an evenly spaced 
% (over time), 2D matrix (size: # of timesteps * # of electrodes), that
% specifies the voltages applied to electrodes over time.
%
% dknapp, 15.8.2023
%
% INPUTS:
%       voltages    :   2D matrix specifying electrode voltages, as above
%       t           :   The time for which interpolation should occur
%       step_times  :   1D matrix with the time for time step.

    start_time = step_times(1);
    end_time   = step_times(end);
    time_steps = length(step_times);
    
    % THIS ASSUMES EQUAL SPACING
    prev_point_idx = ...
        1 + floor((time_steps - 1) * (t - start_time) / (end_time - start_time));
    next_point_idx = prev_point_idx + 1;
    
    interpolated_voltages = ...
        (t - step_times(prev_point_idx)) ...
        * (voltages(next_point_idx, :) - voltages(prev_point_idx, :)) ...
        / (step_times(next_point_idx) - step_times(prev_point_idx)) ...
        + voltages(prev_point_idx, :);
end