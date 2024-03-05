%% Rigid Bodies
% This example follows broadly the same structure as the 'SimpleExample.m'
% script. We additionally create a rigid, rod-shaped collection of ions and
% add them to the simulation.
close all
clearvars

sim = LAMMPSSimulation();
sim.SetSimulationDomain(1e-3,1e-3,1e-3);

% Add a new atom type:
charge = 1;
mass = 40;
Ca = sim.AddAtomType(charge, mass);

% create a cloud of free ions
createIonCloud(sim, 1e-4, Ca, 30);

% position the ions that will form the rigid rod.
rodz = (-2:0.5:2) * 5e-6 + 10e-6; rody = zeros(size(rodz)); rodx = zeros(size(rodz));
rodAtoms = placeAtoms(sim, Ca, rodx', rody', rodz');

% We create a new group, 'rod', which is composed of the rod Atoms, and
% set the Rigid property to true.
rod = sim.Group(rodAtoms);
rod.Rigid = true;

% Create the Paul trap field.
RF = 3.85e6;
z0 = 5.5e-3 / 2;
r0 = 7e-3 / 2;
geometricKappa = 0.244;
U0 = 8;
V0 = 500;
a = -4 * geometricKappa / z0 .^2 / (RF * 2 * pi).^2 * U0 * ((charge * 1.6e-19) / (mass * 1.67e-27));
q = 2 / r0 .^2 / (RF * 2 * pi).^2 * V0 * ((charge * 1.6e-19) / (mass * 1.67e-27));
fprintf('System has a=%.5f, q=%.5f\n', a,q)
sim.Add(linearPT(V0, U0, z0, r0, geometricKappa, RF));

% sim.Add(thermalVelocities(1, 'yes'));
sim.Add(langevinBath(0, 1e-5));

% Output data
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 2));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));
sim.Add(evolve(5000));
sim.Execute();

%% Load data and plot
% Loads simulated trajectories. Plots a small animation showing the rod and
% ions moving together.

% Load results from the output file:
[timestep, id, x,y,z] = readDump('positions.txt');
x = x*1e6; y = y*1e6; z = z*1e6; 

totalIons = size(x, 1);
rodMask = ismember(1:totalIons,[rodAtoms.ID]);

% Determine colors of ions
pastelBlue = [112 146 190]/255;
pastelRed = [237 28 36]/255;
c = repmat(pastelBlue, totalIons, 1);
c(rodMask,:) = repmat(pastelRed, length(rodAtoms), 1);

f = clf;

% Plot animation
while (ishandle(f))
    for i=1:5:size(x,2)
        
        if ~ishandle(f)
            % stop if window is closed.
            break;
        end
        
        depthPlot(x(:,i),y(:,i),z(:,i), c, [100 250]);
        axis equal;
        
        % set axis labels and limits
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
        zlabel('Z (\mum)');
        xlim( [ min(x(:)) max(x(:)) ] );
        ylim( [ min(y(:)) max(y(:)) ] );
        zlim( [ min(z(:)) max(z(:)) ] );
        title(sprintf('Frame %d', i));
        pause(0.04);
    end
end