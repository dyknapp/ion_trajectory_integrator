m = 2.0; % amu
q = 1.0; % electron charges
omega      = 2 * pi * 1.0e+7; % rad/s
trap_potential = 100.0; % V
R          = 5.0; % mm
% TODO: proper Poisson
coll_rate = 1.0e-2; % MHz
room_temp = 300.; % K

max_t      = 1000.0; % us
max_dist   = 1.0e-5; % mm
averaging_window = 100; % samples
pts              = 10000; % samples


position = [1. 1. 1.];
velocity = [0. 0. 0.];

tic
[four_trajectory, its, rpts, collisions] ...
    = inf_trap_bg_gas(position, velocity, m, q, coll_rate, room_temp, ...
                           omega, trap_potential, R, max_t, max_dist, ...
                                              averaging_window, pts);
toc
fprintf("Iterations: %d\n", its);
fprintf("Collisions: %d\n", collisions);

four_trajectory = four_trajectory(1:rpts, :);

plot(four_trajectory(:, 4), four_trajectory(:, 1))