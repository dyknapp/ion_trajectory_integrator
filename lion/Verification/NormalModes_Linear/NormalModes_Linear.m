%% Normal Modes (Linear)
% Simulate a linear string of ions to determine the normal modes. Compare
% with analytical predictions (see: James, App Phys B 66 (2) 1998)
% Author: Elliot Bentine 2019
clear all
close all

% Trap parameters are defined below.
Trapfrequency = 3.85e6; %Hz
TrapQ = 0.3;
TrapA = -0.001;
R0 = (7e-3)/2; %m
EndcapZ0 = (5.5e-3)/2; %m
geometricC = 0.244;

% Simulation results will be stored in the result and output folders.
f = java.io.File('output');
if ~f.exists()
    mkdir('output');
end

f = java.io.File('result');
if ~f.exists()
    mkdir('result');
end

ionMass = 40; %amu
ionCharge = +1; % electron charges
q_z = 0;
a_z = -2*TrapA;
theoreticalZFreq = Trapfrequency/2*sqrt((q_z^2)/2 + a_z);
axialAngularFreq = 2*pi*theoreticalZFreq;
lengthScale = (((ionCharge * 1.6e-19) .^2)/(4*pi*8.85e-12)/(ionMass * 1.661e-27 *(axialAngularFreq.^2))).^(1/3);

NumberOfIons = 1:7;
results = struct('N', {}, 'z', {}, 'time', {});

for i=1:length(NumberOfIons)
    
    % We create a simulation in LAMMPS using the Matlab Wrapper. We have
    % used a stronger damping than would typically be measured in
    % experiment for the Langevin bath in order that the series of
    % simulations may be completed in reasonable time.
    sim = LAMMPSSimulation();
    
    SetSimulationDomain(sim, 1e-3,1e-3,1e-3);
    ions = sim.AddAtomType(ionCharge, ionMass);
    [oscV, endcapV] = getVs4aq(ions, Trapfrequency, EndcapZ0, R0, geometricC, TrapA, TrapQ);
    createIonCloud(sim, 1e-3, ions, NumberOfIons(i), 1e-4);
    sim.Add(linearPT(oscV, endcapV, EndcapZ0, R0, geometricC, Trapfrequency));
    bath = langevinBath(3e-4, 2e-5);
    sim.Add(bath);
    
    % Run simulation to minimise the ions
    sim.Add(evolve(40000));
    
    sim.Remove(bath);
    % Deliver periodic 'kicks' by randomizing the velocity, such that we
    % excite many modes over the simulation. This is performed by the
    % langevin bath.
    sim.Add(langevinBath(3e-4, 1e-4));
    sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 10));
    sim.Add(evolve(100000));
    sim.Execute();
    
    % Load the trajectories for analysis
    [timesteps, ~, x,y,z] = readDump('positions.txt');
    
    x=x';
    y=y';
    z=z';
    
    time = timesteps*sim.TimeStep;
    
    % Plot a graph showing vertical positions, so that we can verify the
    % crystal remains intact.
    clf; plot3(x,y,z); view([0, 0]); axis equal; pause(0.1);
    results(end+1) = struct('N', NumberOfIons(i), 'z', z, 'time', time);
end

% Get spectrum of each normal mode
output = {};
for i=1:length(results)
    result = results(i);
    z = result.z;
    % sort the ions
    [~,indices] = sort(z(1,:));
    z = z(:,indices);
    
    [freqs, modes] = getNormalModeSpectrum(result.N);
    eZ = z * modes';
    eZ = eZ - mean(eZ,1); % remove dc component from each mode
    
    % enumerate through modes
    amplitudes = {};
    for j=1:result.N
        [f,a] = getFT(result.time, eZ(:,j));
        amplitudes{j} = a.*conj(a);
    end
    
    result.amps = cell2mat(amplitudes);
    result.freqs = freqs;
    result.f = f;
    output{end+1} = result;
end
output = [output{:}];


save('results.mat', 'results', 'axialAngularFreq', 'output');