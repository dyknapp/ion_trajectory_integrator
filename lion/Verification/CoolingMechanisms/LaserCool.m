%% Verify Laser Cooling
% Author: Elliot Bentine 2016
% 
% Verify that laser cooling simulations in lion reduces the velocity
% exponentially along the chosen direction but does not alter other
% velocity components.

%% Test 1: Zero-temperature bath, damping coefficients
% Test that the langevin bath damps velocity at the correct rate.

sim = LAMMPSSimulation();
sim.SetSimulationDomain(1e-3,1e-3,1e-3);

% Add a new atom type:
charge = 0;
mass = 40;
ion = sim.AddAtomType(charge, mass);

% Initialise non interacting atoms with a random velocity. Log the velocity
% output over time.

% Chosen random direction to perform cooling along
timeConst = 1e-6;
randDir = rand(3,1);
randDir = randDir / norm(randDir);
gamma = randDir / timeConst;

createIonCloud(sim, 1e-6, ion, 30, 0);
Add(sim, thermalVelocities(1, 'no', 1));
Add(sim, laserCool(ion, gamma));

sim.TimeStep = 1e-9;
Add(sim, dump('v.txt', {'id', 'vx','vy','vz'}, 10));

sim.Add(evolve(1000));
sim.Execute();

%%
% Read data, fit velocity damping and compare to expected value.
[timestep, ~, vx, vy, vz] = readDump('v.txt');

% project velocities along different axes
a = cross(randDir, [ 1 0 0]');
b = cross(randDir, a);
c = randDir;

v = vx .* randDir(1) + vy .* randDir(2) + vz .* randDir(3);
va = vx .* a(1) + vy .* a(2) + vz .* a(3);
vb = vx .* b(1) + vy .* b(2) + vz .* b(3);

% Perform exponential fit for each atom trajectory to determine damping
% constants
t = timestep*sim.TimeStep';
csim = [];
for i=1:size(v, 1)
    f = fit(t', v(i,:)', 'exp1');
    csim(i) = f.b;
end

fprintf('\n\n');
fprintf('Laser cooling test simulation results:\n');
fprintf('--------------------------------------\n');
fprintf('Exponential fit for N=%d atoms:\n', size(v,1))
fprintf('c: mean = %.4e, std = %.4e\n', mean(csim), std(csim))
fprintf('Original value = %.4e\n', -norm(c))

h1 = plot(t*1e6, va', 'Color', [ 0.8 0.8 0.8 ]); hold on
h2 = plot(t*1e6, vb', 'Color', [ 0.5 0.5 0.5 ]); hold on
h3 = plot(t*1e6, v', 'Color' , [ 0.8 0.4 0.4 ]); hold off
xlabel('Time ($\mu$s)', 'Interpreter', 'Latex');
ylabel('Velocity (m/s)');
title('Velocity parallel and perpendicular to laser beam', 'Interpreter', 'Latex');
legend([ h1(1), h2(1), h3(1) ], { 'perp 1', 'perp 2', 'parallel' });
axis tight