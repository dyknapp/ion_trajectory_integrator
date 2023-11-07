% Add r = -1 and r = -2 rows to the potential arrays so that we can still
% treat r = 0 normally without worrying about going out of bounds.
% Unfortunately, in MATLAB we aren't able to do this with the native array
% syntax.  It works more naturally in FORTRAN.  For MATLAB, r = -2 will be
% element 1 vice versa.
%
% dknapp, 01.11.2023
%
%
% INPUT:
%       


function output = prep_cylindrical(potential_maps)
    dimensions = size(potential_maps);
    output = zeros(dimensions + [0 2 0]);
    output(:, 3:end, :) = potential_maps;
    output(:, 2, :) = potential_maps(:, 2, :); % r = +/- 1
    output(:, 1, :) = potential_maps(:, 3, :); % r = +/- 2
end