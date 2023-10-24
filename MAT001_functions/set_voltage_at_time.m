% Function for manipulating the time-dependent voltage arrays used as input
% for the trajectory calculation.
% At the specified time, the electrode with the specified index will be
% switched to the specified voltage.
%
% See MAT007_pure_matlab_examples/electrode_stack_cartesian.m
%
% 24.10.2023   dknapp: wrote the function

function voltages_array = set_voltage_at_time(electrode, voltage, time, step_times, voltages_array)
    t_start_idx = find(time <= step_times, 1);
    voltages_array(t_start_idx:end, electrode) = voltage;
end