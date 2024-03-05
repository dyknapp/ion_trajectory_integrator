%% Verify Langevin bath
% Author: Elliot Bentine 2016
% 
% Verify that the langevin bath in LIon gives an isotropic exponential
% damping to velocity for zero temperature bath. Verify that for a finite
% temperature bath the ensemble thermalises.
% 
% Create N random ions, all at initialT. Each is coupled to a separate
% Langevin bath, with a different time constant and temperature.
% 
% The simulation is performed simultaneously for multiple sets of N ions
% each.

initialT = 15 * 1e-6;
bathT = [ 1 5 10 20 30 40 ] * 1e-6;
timeConstant = [ 5 15 10 20 10 5 ] * 1e-6;

sim = LAMMPSSimulation();
sim.TimeStep = 1e-8;
sim.SetSimulationDomain(1e-3,1e-3,1e-3);

charge = 0;
mass = 40;
ionType = sim.AddAtomType(charge, mass);

% Add groups of N ions to the simulation
N = 100;
groups = {};
ions = {};
for i=1:length(bathT)
    x = (1:N)*1e-6;
    y = x;
    z = x + i*1e-6;
    ions{i} = placeAtoms(sim, ionType, x', y', z');
    groups{i} = sim.Group(ions{i});
end

startBath = langevinBath(initialT, 1e-6);
sim.Add(startBath);
sim.Add(evolve(100));
sim.Remove(startBath);
sim.Add(thermalVelocities(initialT, 'no', 100));

% Add a bath for each group of ions
for i=1:length(ions)
    sim.Add(langevinBath(bathT(i), timeConstant(i), groups{i}, i)); 
end

Add(sim, dump('v.txt', {'id', 'vx','vy','vz'}, 10));

sim.Add(evolve(5000));
sim.Execute();

% Read the velocities, determine temperature of each group.
[timestep, ~, vx, vy, vz] = readDump('v.txt');
ke = (vx.^2 + vy.^2 + vz.^2) * mass * Const.amu / 2;
groupKe = cellfun(@(ionGroup) sum(ke([ionGroup.ID], :), 1), ions, 'UniformOutput', 0);
groupKe = cat(1, groupKe{:});
T = 2/3 * movmean(groupKe', 3) / N / Const.kB;
t = sim.TimeStep * timestep;
t = t-min(t);

save('results.mat', 't', 'T', 'initialT', 'bathT', 'timeConstant');
clf; 
plot(t,T); hold on;
% Plot theoretical lines for each group. It should be a decay from initial
% temperature to final temperature with rate constant 2*c (because T
% \propto v^2).
for i=1:length(bathT)
    plot(t, initialT + (bathT(i) - initialT)*(1-exp(1).^(-2*t/timeConstant(i))), 'k:');
end