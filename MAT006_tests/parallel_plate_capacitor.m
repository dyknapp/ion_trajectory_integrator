%% Variables to set beforehand
% For this example, we will just use an analytical expression for the
% potentials

% Time settings in microseconds
start_time =  0.0;
end_time   = 25.0;

% Number of particles
N_ions = 1;

% velocity verlet settings
m = 1.0;                % amu
q = 1.0;                % atomic units
d = 1.0;                % mm
maxdist =  0.01;      % mm

time_steps = 1000;
step_times = linspace(start_time, end_time + ((end_time - start_time) / 100), time_steps);
voltages = zeros([time_steps, 2]);

%% RF and endcap electrode voltages.
% We just have a parallel place capacitor, using the potentials defined in
% the next section of the code.  Our two electrodes will be at the z
% extremes.  Electrode 1 is at z = 0;

voltages(:, 1) =  50.3 * ones([time_steps 1]);
voltages(:, 2) = -50.3 * ones([time_steps 1]);

%% Loading & Reading
dimensions = [51, 51, 1006];

% According to Laplace's equation, the potential needs to be linear in
% vacuum so, setting electrode 1 to 1V and electrode 2 to 0V, we get,
potential_maps = zeros([2 dimensions]);
for z = double(1:dimensions(3))
    potential_maps(1, :, :, z) = (1.0 - z / double(dimensions(3))) * ones(dimensions(1:2));
end

% And the same the other way around
for z = double(1:dimensions(3))
    potential_maps(2, :, :, z) = z / double(dimensions(3)) * ones(dimensions(1:2));
end

% Invisible electrodes
is_electrode = zeros(dimensions);

%% Initializing

xx1 = 10.0;
yy1 = 10.0;
zz1 = 10.0;

vxx1 = 0.1;
vyy1 = 0.2;
vzz1 = 0.0;

%% Integration: MATLAB version
tic
ion_trajectory =  integrate_trajectory(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                                      potential_maps, voltages, step_times, ...
                                      dimensions, is_electrode, m, q, d, ...
                                      maxdist, end_time);
x_traj = ion_trajectory.x;
y_traj = ion_trajectory.y;
z_traj = ion_trajectory.z;
ts     = ion_trajectory.t;
exs    = ion_trajectory.ex;
eys    = ion_trajectory.ey;
ezs    = ion_trajectory.ez;
elapsed_time = toc;
fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(length(ts) / elapsed_time));

%% Integration: FORTRAN version
if check_FORTRAN_integration_args(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                                  potential_maps, voltages, step_times, ...
                                  time_steps, dimensions, int32(is_electrode), ...
                                  2, m, q, d, maxdist, end_time)
    tic
    potential_maps_size = size(potential_maps);
    potential_maps = reshape(potential_maps, [potential_maps_size(1), dimensions]);
    [x_traj, y_traj, z_traj, ts, exs, eys, ezs, its] ...
        = trajectory_integration_module(xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                          potential_maps, voltages, step_times, ...
                          time_steps, dimensions, int32(is_electrode), ...
                          2.0, m, q, d, maxdist, end_time);
    
    its = int32(its);
    original_length = length(ts);
    x_traj = x_traj(1:its);
    y_traj = y_traj(1:its);
    z_traj = z_traj(1:its);
    ts     =     ts(1:its);
    exs    =    exs(1:its);
    eys    =    eys(1:its);
    ezs    =    ezs(1:its);
    
    elapsed_time = toc;
    fprintf("Simulation took %.3gs (%d it/s)\n", elapsed_time, round(its / elapsed_time));

else
    fprintf("Unsafe input detected.\n")
end










