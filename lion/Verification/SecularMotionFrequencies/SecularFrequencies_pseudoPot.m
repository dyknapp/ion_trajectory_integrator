%% Verify Secular Motion Frequency ( pseudopotential )
% Author: Elliot Bentine 2013
% 
% This script verifies that the linear Paul trap pseudopotential
% implementation in (py)Lion gives the correct force. It uses a single ion
% and is as such very fast to run, but as a result it cannot verify that
% ion-ion interactions are correct. Mathieu parameters for the
% translational motion (a and q) are varied and the secular frequencies
% extracted from the simulation results. These are plotted against
% theoretical predictions.

clear
close all

%define the ions to use.
ionMass = 40;
ionCharge = +1;

%Parameters for the trap.
TrapA = -0.001;
TrapQs = 0.05:0.025:0.3;

%Physical parameters of the trap.
Trapfrequency = 10e6; %Hz
R0 = 1e-3; %m
EndcapZ0 = 1e-3; %m
geometricC = 0.5;

%Arrays to store output values.
pseudoX = zeros(0);
pseudoY = zeros(0);
pseudoZ = zeros(0);

%Iterate through different values of trap q.
for TrapQ = TrapQs

    
%Create an experiment in lammps
sim = LAMMPSSimulation();
SetSimulationDomain(sim, 1e-3,1e-3,1e-3);

C40Ion = AddAtomType(sim, ionCharge, ionMass);
radiusofIonCloud = 1e-3;
createIonCloud(sim, radiusofIonCloud, C40Ion, 1, 1e-4);

ionMasskg = ionMass * 1.66e-27; %in kilograms
ionChargeC = ionCharge * 1.6e-19; %in Coulombs

[oscV, endcapV] = getVs4aq(C40Ion, Trapfrequency, EndcapZ0, R0, geometricC, TrapA, TrapQ);
sim.Add(linearPseudoPT(oscV, endcapV, EndcapZ0, R0, geometricC, Trapfrequency, C40Ion));
    
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 1));

% Run simulation
sim.Add(evolve(10000));
sim.Execute();

%Load the output and analyse it for secular frequencies.
[timesteps, ~, x,y,z] = readDump('positions.txt');
time = (timesteps)*sim.TimeStep;

%%
%Find highest power freq component of FFT of each dimension.
[xFreq, xAmp] = getFT(time, x);
[yFreq, yAmp] = getFT(time, y);
[zFreq, zAmp] = getFT(time, z);

xPow = xAmp .* conj(xAmp);
yPow = yAmp .* conj(yAmp);
zPow = zAmp .* conj(zAmp);

[~, ix] = max(xPow);
[~, iy] = max(yPow);
[~, iz] = max(zPow);

%store freq of secular motion
pseudoX(end+1) = xFreq(ix);
pseudoY(end+1) = yFreq(iy);
pseudoZ(end+1) = zFreq(iz);

end

%Calculate theoretical secular frequencies.
theoreticalRadialFreq = Trapfrequency/2*((TrapQs.^2)/2 + -abs(TrapA)).^0.5;
theoreticalZFreq = Trapfrequency/2*sqrt(-2*TrapA)*ones(size(TrapQs));

%Plot measured secular frequencies from LAMMPS as crosses, and theoretical
%frequencies from theory as circles.
gr = figure; hold on;
plot(TrapQs, 1e-3*pseudoX, 'rx');
plot(TrapQs, 1e-3*pseudoY, 'gx');
plot(TrapQs, 1e-3*pseudoZ, 'bx');
plot(TrapQs, 1e-3*theoreticalRadialFreq, 'ro');
plot(TrapQs, 1e-3*theoreticalRadialFreq, 'go');
plot(TrapQs, 1e-3*theoreticalZFreq, 'bo');

title('Pseudopotential Trap Frequencies for a Single Ion');
xlabel('Trap Q');
ylabel('Secular freq (kHz)');
legend('LAMMPS x', 'LAMMPS y', 'LAMMPS z', 'Theory x', 'Theory y', 'Theory z');
saveas(gr, 'pseudopot.fig');
save('results.mat', ...
    'TrapQs', ...
    'pseudoX', 'pseudoY', 'pseudoZ', ...
    'theoreticalRadialFreq', ...
    'theoreticalZFreq' ...
    );