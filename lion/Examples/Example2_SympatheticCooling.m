%% Example: Sympathetic Cooling
% Author: Elliot Bentine.
%
% This simulations follows in the spirit of an earlier paper by Zhang,
% Offenberg, Roth, Wilson, Schiller (PRA 76 012719 2007) which can be found
% at this url:
% http://www.spektron.de/Publikationen/2007/PhysRevA_76_012719.pdf
%
% A two species Coulomb crystal is formed by coupling both species to a 1mK
% Kelvin temperature bath. At a specified time t1, the bath is removed and
% both species heat via micromotion in the trap. At a second time t2,
% cooling is enabled for one species, sympathetically cooling the second.

% Define timesteps
minimisationSteps = 100000;
interval = 60000;

% Define trap parameters
rf = 5.634e6; % Hz
Vo = 252.5; % V
Ve = 8; % V
geomC = 0.3;
r0  = 3.5e-3;
z0  = 10e-3;

% Create the simulation instance
sim = LAMMPSSimulation();
sim.GPUAccel = 0;
SetSimulationDomain(sim, 1e-3,1e-3,1e-3);

% The two ion species used in this paper are NH3+ and 40Ca+
NH3 = AddAtomType(sim, 1, 17);
Ca40 = AddAtomType(sim, 1, 40);

% Create the ion clouds.
radius = 5e-4; % place atoms randomly within this radius
NumberNH3 = 20;
NumberCa  = 50;
NH3Ions = createIonCloud(sim, radius, NH3, NumberNH3);
Ca40Ions = createIonCloud(sim, radius, Ca40, NumberCa);

% add the electric field
sim.Add(linearPT(Vo, Ve, z0, r0, geomC, rf));

% minimise the system using a langevin bath applied to both species
allBath = langevinBath(1e-3, 30e-6);
sim.Add(allBath);
sim.Add(evolve(minimisationSteps));

% Having minimised the system, output the coordinates of both species
% and time averaged velocities from which we can calculate the secular
% temperature.
sim.Add(dump('sympcool.txt', {'id', 'x', 'y', 'z', timeAvg({'vx', 'vy', 'vz'}, 1/rf)}, 20));
sim.Add(evolve(interval));

% At t1, remove the langevin bath and add laser cooling for just one species
sim.Remove(allBath);
sim.Add(laserCool(Ca40, [ 1e5 0 0 ]));
sim.Add(evolve(interval*2));

sim.Execute();

%% Load results
% Load simulation results and calculate the temperature of each
% species.

[t, id, x,y,z, sx,sy,sz] = readDump('sympcool.txt');
t = (t-minimisationSteps)*sim.TimeStep;

v2 = @(ind) sum(sx(ind, :).^2 + sy(ind, :).^2 + sz(ind,:).^2, 1);
T_NH3 = v2([NH3Ions.ID]) * Const.amu * NH3.Mass / 3 / Const.kB;
T_Ca = v2([Ca40Ions.ID]) * Const.amu * Ca40.Mass / 3 / Const.kB;

color_NH3 = [112 146 190]/255;
color_Ca = [237 28 36]/300;255;

clf;

% Plot a representation of the final Coulomb crystal
ax = subplot(1, 4, 1);
color = repmat(color_Ca, size(x, 1), 1);
color([NH3Ions.ID], :) = repmat(color_NH3, length([NH3Ions.ID]), 1);
depthPlot(1e6*x(:,end),1e6*y(:,end),1e6*z(:,end),color, [ 40 80 ]);
axis equal
view(45, 22);

% Plot the energies of the two species as a function of time
subplot(1, 4, [ 2 3 4 ]);
plot(t*1e6, T_NH3 / NumberNH3 * 1e3, '-', 'Color', color_NH3*0.8); hold on;
plot(t*1e6, T_Ca / NumberCa * 1e3, '-', 'Color', color_Ca); hold off
xlabel('time (\mus)');
ylabel('Temperature (mK)');
