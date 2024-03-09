clearvars
constants = physical_constants();

% Simulation timescales (us)
end_time   =  100.0;
record_step =  0.001;
burst_time =   0.000;

% Simulation step scale (mm)
maxdist = 1.0e-8;

%% Initial conditions
orbit_omega = 2 * pi * 1.0e+6; % rad/s
d = (constants("elementary charge")^2 / (2 * pi * constants("epsilon0") * constants("proton mass") * orbit_omega^2))^(1/3);

particles = 2;
% mm, mm/us
positions  = (1.0e+3) * [ 0 0; [d/2. -d/2.]; 0 0]';
velocities = (1.0e-3) * [orbit_omega*d/2 -orbit_omega*d/2; 0 0; 0 0]';
ms = [1  1];
qs = [1 -1];

% No ion trap potential
omega = 1.0;
depth = 0.0;
R = 1.0;

%% Simulation

mh = mexhost;
tic;
[trajectories, times, its, recorded] = ...
    feval(mh, 'nbody', ...
        particles, ...
        positions, ...
        velocities,  ...
        ms,  ...
        qs, ...
        omega,  ...
        depth,  ...
        R,  ...
        end_time, ...
        maxdist, ...
        record_step, ...
        burst_time ...
    );
elapsed = toc;
fprintf("Simulation finished ( %5.1f s).  Iterations: %9.1f ( %.3g its/s). Mean timestep: %.3g us, Recorded points: %d\n", elapsed, its, its / elapsed,  (times(recorded) / double(its)), recorded)

% %%
% plot(squeeze(trajectories(1:recorded, 2, 1)), squeeze(trajectories(1:recorded, 2, 2)), '.-')
% axis image

%%

error = sqrt((squeeze(trajectories(1:recorded, 1, 1)) - (sin(orbit_omega*times*(1.0e-6))*(1.0e+3)*d/2.)).^2 ...
           + (squeeze(trajectories(1:recorded, 1, 2)) - (cos(orbit_omega*times*(1.0e-6))*(1.0e+3)*d/2.)).^2);
loglog(times, error)
set(gca,'TickDir','out');  grid("on");