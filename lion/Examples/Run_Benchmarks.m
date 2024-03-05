%%
% Runs the benchmark simulation for different atom numbers and records the
% wall time required to complete the simulation.

Numbers = [ 1 2 3 6 10 18 32 56 100 178 316 562 800 1000 2000 ];
GPUAccel = 0;
results = struct('N', {}, 'time', {});

for NumberOfIons=Numbers
    wallTime = Benchmark(NumberOfIons, GPUAccel);
    results(end+1) = struct('N', NumberOfIons, 'time', wallTime);
end

save('bench.mat', 'Numbers', 'GPUAccel', 'results');