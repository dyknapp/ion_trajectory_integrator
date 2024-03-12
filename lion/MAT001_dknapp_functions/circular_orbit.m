clearvars
constants = physical_constants();

% Simulation timescales (us)
end_time   =   100000.0;
record_step =  2.0;
burst_time =   0.00;

% Simulation step scale (mm)
maxdist = 1.0e-3;

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
fprintf("Simulation finished (VV) ( %5.1f s).  Iterations: %9.1f ( %.3g its/s). Mean timestep: %.3g us, Recorded points: %d\n", elapsed, its, double(its) / elapsed,  (times(recorded) / double(its)), recorded)
tic;
[trajectories4, times4, its4, recorded4] = ...
    feval(mh, 'nbody4', ...
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
% [trajectories, times, its, recorded] = ...
%   nbody(...
%         particles, ...
%         positions, ...
%         velocities,  ...
%         ms,  ...
%         qs, ...
%         omega,  ...
%         depth,  ...
%         R,  ...
%         end_time, ...
%         maxdist, ...
%         record_step, ...
%         burst_time ...
%     );
elapsed = toc;
fprintf("Simulation finished (Y4) ( %5.1f s).  Iterations: %9.1f ( %.3g its/s). Mean timestep: %.3g us, Recorded points: %d\n", elapsed, its, double(its) / elapsed,  (times(recorded) / double(its)), recorded)

% %%
% plot(squeeze(trajectories(1:recorded, 2, 1)), squeeze(trajectories(1:recorded, 2, 2)), '.-')
% hold on
% plot(squeeze(trajectories(1:recorded, 1, 1)), squeeze(trajectories(1:recorded, 1, 2)), '.-')
% hold off
% axis image

%%

error = sqrt((squeeze(trajectories(1:recorded, 1, 1)) - (sin(orbit_omega*times(1:recorded)*(1.0e-6))*(1.0e+3)*d/2.)).^2 ...
           + (squeeze(trajectories(1:recorded, 1, 2)) - (cos(orbit_omega*times(1:recorded)*(1.0e-6))*(1.0e+3)*d/2.)).^2) / ((1.0e+3)*d/2.);
error4 = sqrt((squeeze(trajectories4(1:recorded, 1, 1)) - (sin(orbit_omega*times4(1:recorded)*(1.0e-6))*(1.0e+3)*d/2.)).^2 ...
            + (squeeze(trajectories4(1:recorded, 1, 2)) - (cos(orbit_omega*times4(1:recorded)*(1.0e-6))*(1.0e+3)*d/2.)).^2) / ((1.0e+3)*d/2.);
loglog(times(1:recorded), error, '.-')
hold on
loglog(times(1:recorded), error4, '.-')
hold off
set(gca,'TickDir','out');  grid("on");
legend("Velocity Verlet", "Yoshida 4")

xlabel('Time (us)')
ylabel('relative position error')
title('1MHz circle Long Time Position Error Comparison (md=1e-3)')

%% 
k = (constants("elementary charge")^2) / (4 * pi * constants("epsilon0"));
distances  = vecnorm(squeeze( trajectories(1:recorded, 1, 1:3)) -  squeeze(trajectories(1:recorded, 2, 1:3)), 2, 2);
distances4 = vecnorm(squeeze(trajectories4(1:recorded, 1, 1:3)) - squeeze(trajectories4(1:recorded, 2, 1:3)), 2, 2);
energies  = ((0.5 * constants("proton mass")) * vecnorm(squeeze(trajectories(1:recorded, 1, 4:6)), 2, 2).^2)  + (k / distances)';
energies4 = ((0.5 * constants("proton mass")) * vecnorm(squeeze(trajectories4(1:recorded, 1, 4:6)), 2, 2).^2) + (k / distances4)';
E0 = (0.5 * constants("proton mass")) * norm(velocities(1, :))^2.;
loglog(times(1:recorded), abs(energies - E0) / E0)
hold on
loglog(times4(1:recorded), abs(energies4 - E0) / E0)
hold off
set(gca,'TickDir','out');  grid("on");
legend("Velocity Verlet", "Yoshida 4")

xlabel('Time (us)')
ylabel('relative energy error')
title('1MHz circle long time energy error comparison (md=1e-3)')

%%

loglog(times(1:recorded), vecnorm(squeeze(trajectories(1:recorded, 1, 4:6)) + squeeze(trajectories(1:recorded, 2, 4:6)), 2, 2))
hold on
loglog(times4(1:recorded), vecnorm(squeeze(trajectories4(1:recorded, 1, 4:6)) + squeeze(trajectories4(1:recorded, 2, 4:6)), 2, 2))
hold off
legend("Velocity Verlet", "Yoshida 4")
xlabel('Time (us)')
ylabel('momentum error')
title('1MHz circle long time momentum error comparison (md=1e-3)')