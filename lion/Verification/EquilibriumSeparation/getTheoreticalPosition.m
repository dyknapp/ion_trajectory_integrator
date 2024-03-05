function [ positions , success ] = getTheoreticalPosition( numberOfAtoms, lengthScale )
%GETTHEORETICALPOSITION Gets theoretical positions of ions in a linear paul
%trap when a chain configuration is formed. Adapted from 'Quantum dynamics
%of cold trapped ions with application to quantum computation' by DFV
%James, App. Phys. B 1998

positions = [];

switch numberOfAtoms
    case 1
        positions = 0;
    case 2
        positions = [-0.62996 0.62996];
    case 3
        positions = [-1.0772 0 1.0772];
    case 4
        positions = [-1.4368 -0.45438 0.45438 1.4368];
    case 5
        positions = [-1.7429 -0.8221 0 0.8221 1.7429];
    case 6
        positions = [ -2.0123 -1.1361 -0.36992 0.36992 1.1361 2.0123];
    case 7
        positions = [-2.2545 -1.4129 -0.68694 0 0.68694 1.4129 2.2545];
    case 8
        positions = [-2.4758 -1.6621 -0.96701 -0.31802 0.31802 0.96701 1.6621 2.4758];
    case 9
        positions = [-2.6803 -1.8897 -1.2195 -0.59958 0 0.59958 1.2195 1.8897 2.6803];
    case 10
        positions = [-2.8708 -2.10003 -1.4504 -0.85378 -0.2821 0.2821 0.85378 1.4504 2.10003 2.8708];     
end

success = ~isempty(positions);

positions = positions*lengthScale;


end