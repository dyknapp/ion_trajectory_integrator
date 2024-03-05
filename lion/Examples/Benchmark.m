function walltime = Benchmark(NumberOfIons, GPUAccel)

if nargin < 1
    NumberOfIons = 1000;
end

if nargin < 2
    GPUAccel = 1;
end

sim = LAMMPSSimulation();
sim.GPUAccel = GPUAccel;
sim.Verbose = 0;

sim.SetSimulationDomain(1e-3,1e-3,1e-3);

% Add some atoms:
charge = 1;
mass = 100;
ions = sim.AddAtomType(charge, mass);
createIonCloud(sim, 1e-5, ions, NumberOfIons, 1337);

% Add trap
RF = 3.85e7;
z0 = 5.5e-3; r0 = 7e-3; geomC = 0.244;
[vs(1), vs(2)] = getVs4aq(ions, RF, z0, r0, geomC, -0.01, 0.3);
sim.Add(linearPT( vs(1), vs(2), 5.5e-3, 7e-3, 0.244, RF));

% Add bath
sim.Add(langevinBath(3e-3, 1e-5));

%Configure outputs.
sim.Add(dump('positions.txt', {'id', 'x', 'y', 'z'}, 100));

% Run simulation
sim.Add(evolve(100000));

timer = tic;
sim.Execute();
walltime = toc(timer);

figure(1); clf; set(gcf, 'Color', 'w');
[~, ~, x,y,z] = readDump('positions.txt');
plot3(x',y',z','-','Color', [0.8 0.8 0.8]); hold on
plot3(x(:,end)',y(:,end)',z(:,end)','k.');
plot3(0,0,0,'g+');  hold off
end